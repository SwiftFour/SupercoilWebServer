/* Author: Qin Ma <maqin@csbl.bmb.uga.edu>, Jan. 25, 2010
 * Usage: This is part of the supercoil package. Use, redistribute, modify
 * without limitations.
 *
 * Produces two graphs sequentially, derived from microarray data.
 */

#include "expression_analysis.h"
#include <vector>
#include <algorithm>
 /*we can reduce the HEAP_SIZE when the data contain so many genes so that memory is not enough*/
static const int HEAP_SIZE = 20000000;

/**************************************************************************/
/* String intersection function without string copying, only numbers */
/*caculate the weight of the edge in the first graph*/
static int str_intersect_r(const discrete *s1, const discrete *s2)
{
  int common_cnt = 0;
  /* s1 and s2 of equal length, so we check s1 only */
  int i;
  for (i = 0;i < cols;i++)
  {
    if (*s1 == *s2 && (*s1 != 0))
      common_cnt++;
    s1++;
    s2++;
  }
  return common_cnt;
}

/**************************************************************************/
continuous get_spearman(const discrete *s1, const discrete *s2, int row_1, int row_2, int cnt)
{
  std::vector<discrete> ss1(cnt), ss2(cnt);
  std::vector<continuous> ss11(cnt), ss22(cnt), temp1(cnt), temp2(cnt);
  continuous score = 0, ave1 = 0, ave2 = 0, var1 = 0, var2 = 0;
  int i = 0, j = 0;
  for (i = 0;i < cols;i++)
  {
    if (s1[i] == s2[i] && (s1[i] != 0))
    {
      ss11[j] = arr[row_1][i];
      ss22[j] = arr[row_2][i];
      temp1[j] = arr[row_1][i];
      temp2[j] = arr[row_2][i];
      j++;
    }
    i++;
  }
  std::sort(temp1.begin(), temp1.end());
  for (i = 0;i < cnt;i++)
  {
    ss1[i] = 0;
    for (j = 0;j < cnt;j++)
      if (ss11[i] == temp1[j])
        ss1[i] = j;
  }
  std::sort(temp2.begin(), temp2.end());
  for (i = 0;i < cnt;i++)
  {
    ss2[i] = 0;
    for (j = 0;j < cnt;j++)
      if (ss22[i] == temp2[j])
        ss2[i] = j;
  }
  /*get var and ave*/
  for (j = 0; j < cnt; j++)
  {
    ave1 += ss1[j];
    ave2 += ss2[j];
  }
  ave1 = ave1 / cnt;
  ave2 = ave2 / cnt;
  for (j = 0; j < cnt; j++)
  {
    var1 += (ss1[j] - ave1)*(ss1[j] - ave1);
    var2 += (ss2[j] - ave2)*(ss2[j] - ave2);
  }
  var1 = sqrt(var1);
  var2 = sqrt(var2);
  for (i = 0; i < cnt; i++)
    score += (ss1[i] - ave1)*(ss2[i] - ave2);
  score = fabs(score / (var1*var2));
  return score;
}
/*************************************************************/
/* String intersection function without string copying, only numbers */
/*caculate the weight of the edge with negative regulation in the first graph*/
static int str_negative_intersect_r(const discrete *s1, const discrete *s2)
{
  int common_cnt = 0;
  /* s1 and s2 of equal length, so we check s1 only */
  int i;
  for (i = 0;i < cols;i++)
  {
    if ((s1[i] != 0) && (symbols[s1[i]] == -symbols[s2[i]]))
      common_cnt++;
  }
  return common_cnt;
}
/**************************************************************************/
/*track the gene id in ptt file*/
int which_gene_in_ptt(int num)
{
  int i, id = 0;
  for (i = 0; i < geneNum; i++)
  {
    if (strcmp(genes_n[num], gene[i]->synonym) == 0)
    {
      id = i + 1;
      break;
    }
  }
  return id;
}

/**************************************************************************/

void expression_parse(const char* fn, const char* fn1, const char* fn2, const char* fn3, const char* fn4)
{
  FILE *fw = mustOpen(fn, "w");
  FILE *fw1 = mustOpen(fn1, "w");
  FILE *fw2 = mustOpen(fn2, "w");
  FILE *fw3 = mustOpen(fn3, "w");
  FILE *fw4 = mustOpen(fn4, "w");
  int i;
  int j;
  int k;
  int cnt;
  int cnt1;
  int cnt2;
  int g1;
  int g2;
  double distance = 0;
  if (po->COL_WIDTH == 2)
    po->COL_WIDTH = MAX(cols / 20, 2);

  /*the virable for spearman calculate*/
  continuous spearman = 0; int cnt_r = 0;

  /* Generating seed list and push into heap */
  /*progress("Generating co-expression gene pairs (minimum weight %d)", po->COL_WIDTH); */
  std::vector<int> cache(rows);
  for (i = 0; i < rows; i++)
    cache[i] = which_gene_in_ptt(i);
  /* iterate over all genes to retrieve all qualified co-expression gene pairs and get lost function*/
  for (i = 0; i < rows; i++)
  {
    if (i % 10 == 0)
      verboseDot();
    for (j = i + 1; j < rows; j++)
    {
      cnt = str_intersect_r(arr_c[i], arr_c[j]);
      cnt1 = str_intersect_r(arr_c[i], arr_c[i]);
      cnt2 = str_intersect_r(arr_c[j], arr_c[j]);

      if (cnt > po->COL_WIDTH && (abs(cache[i] - cache[j]) < 20 || abs(cache[i] - cache[j]) > (geneNum - 20)))
      {
        cnt_r = str_negative_intersect_r(arr_c[i], arr_c[j]);
        cnt = MAX(cnt, cnt_r);
        /*get spearman rank corelation*/
        if (cnt > (cnt1*cnt2) / cols)
        {
          spearman = get_spearman(arr_c[i], arr_c[j], i, j, cnt);
          g2 = MIN(cache[i], cache[j]);
          if (g2 > 0 && ((cnt > po->CNT && spearman > po->SPEARMAN) || spearman > po->SPEARMAN + .2))
          {
            g1 = MAX(cache[i], cache[j]);
            if (g1 - g2 <= 20)
            {
              for (k = g2 - 1; k < g1 - 1; k++)
              {
                distance = MIN(k - g2 + 2, g1 - k - 1);
                expression_score[k] += cnt*spearman / distance;
                expression_score_1[k] += 1 / distance;
                expression_num[k] ++;
              }
            }
            else if (g1 - g2 >= geneNum - 20)
            {
              distance = geneNum + g2 - g1;
              for (k = g1 - 1; k < geneNum; k++)
              {
                expression_score[k] += cnt*spearman / distance;
                expression_score_1[k] += 1 / distance;
                expression_num[k] ++;
              }
              for (k = 0; k < g2 - 1; k++)
              {
                expression_score[k] += cnt*spearman / distance;
                expression_score_1[k] += 1 / distance;
                expression_num[k] ++;
              }
            }
            fprintf(fw1, "%d\t%d\t%d\t%ld\t%ld\t%d\t%ld\t%ld\t%d\t%d\t%.2f\n", i, j, g2, gene[g2 - 1]->start, gene[g2 - 1]->end, g1, gene[g1 - 1]->start, gene[g1 - 1]->end, g1 - g2, cnt, spearman);
          }
        }
      }
    }
  }
  uglyTime("\n4: co-expression details are writen to %s", fn1);

  /*print the caculated lost score in above part*/
  for (i = 0; i < geneNum - 1; i++)
    /*if (gene[i]->IS_operonEnd)*/
    fprintf(fw, "%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%ld\t%ld\t%.1f\t%.2f\t%.2f\t%s\t%s\n",
      i + 1,
      gene[i]->operon_id,
      gene[i]->synonym,
      check_heg(i + 1, 0),
      check_heg(i + 1, 1),
      gene[i + 1]->synonym,
      check_heg(i + 2, 0),
      check_heg(i + 2, 1),
      gene[i]->end,
      gene[i + 1]->start,
      expression_num[i],
      expression_score[i],
      expression_score_1[i],
      check_nap(i + 1, 0),
      check_nap(i + 1, 1));

  uglyTime("5: expression lost function are writen to %s", fn);

  /* smooth the expression values consider left and right point*/
  std::vector<double> expression_score_smooth(operonNum);
  std::vector<double> expression_score_smooth_1(operonNum);
  std::vector<double> expression_score_operon(operonNum);
  std::vector<double> expression_score_operon_1(operonNum);
  j = 0;
  for (i = 0; i < geneNum; i++)
  {
    if (gene[i]->IS_operonEnd)
    {
      expression_score_operon[j] = expression_score[i];
      expression_score_operon_1[j] = expression_score_1[i];
      j++;
    }
  }

  expression_score_smooth[0] = (expression_score_operon[0] + expression_score_operon[1]) / 2;
  expression_score_smooth_1[0] = (expression_score_operon_1[0] + expression_score_operon_1[1]) / 2;
  for (i = 1; i < operonNum - 2; i++)
  {
    expression_score_smooth[i] = (expression_score_operon[i - 1] + expression_score_operon[i] + expression_score_operon[i + 1]) / 3;
    expression_score_smooth_1[i] = (expression_score_operon_1[i - 1] + expression_score_operon_1[i] + expression_score_operon_1[i + 1]) / 3;
  }
  expression_score_smooth[operonNum - 2] = (expression_score_operon[operonNum - 3] + expression_score_operon[operonNum - 2]) / 2;
  expression_score_smooth_1[operonNum - 2] = (expression_score_operon_1[operonNum - 3] + expression_score_operon_1[operonNum - 2]) / 2;


  /*detect the valleys in the lost function*/
  detectValleyG(expression_score, fw2);
  detectValleyG(expression_num, fw2);
  detectValleyG(expression_score_1, fw2);
  uglyTime("6.1: Detected valleys base on expression are writen to %s", fn2);
  detectValleyO(expression_score_operon, fw4);
  uglyTime("6.2: Detected valleys base on expression are writen to %s", fn4);

  /*print the nap information*/
  if (po->IS_add_nap)
  {
    for (i = 0; i < napNum; i++)
      if (nap[i]->IS_interGene)
        fprintf(fw3, "%d\t%s\t%d\t%d\t%ld\n", nap[i]->id, nap[i]->name, nap[i]->gene_id, nap[i]->operon_id, nap[i]->position);
    uglyTime("7: nap binding sites in inter-genetic regions are writen to %s", fn3);
  }
}


/***************************************************************************/
double getDelta(const std::vector<double>& expression)
{
  double min = 100000;
  double max = 0;
  for (auto it = expression.begin(); it != expression.end(); ++it)
  {
    if (*it > max)
      max = *it;
    if (*it < min)
      min = *it;
  }
  return max - min;
}

/* detect valleys on gene level*/
void detectValleyG(const std::vector<double>& expression_score_input, FILE *fw)
{
  double delta;
  int kk = 0;
  int i = 0;
  for (delta = getDelta(expression_score_input) / 10; delta < getDelta(expression_score_input) / 2; delta += 3 * getDelta(expression_score_input) / 100)
  {
    fprintf(fw, "Delta-expression_score: %.2f\n", delta);
    int emi_peaks[MAX_PEAK];
    int absorp_peaks[MAX_PEAK];
    int emi_count = 0;
    int absorp_count = 0;
    int emission_first = 0;
    if (detect_peak(expression_score_input, geneNum - 1, emi_peaks, &emi_count, MAX_PEAK, absorp_peaks, &absorp_count, MAX_PEAK, delta, emission_first))
    {
      printf("There are too many peaks.\n");
      exit(1);
    }

    for (i = 0; i < absorp_count; ++i)
    {
      if (po->IS_add_nap)
        kk = which_nap(absorp_peaks[i] + 1, absorp_peaks[i] + 2, expression_score_input);
      if (kk == 0 || !po->IS_add_nap)
        kk = which_smallest(absorp_peaks[i] + 1, absorp_peaks[i] + 2, expression_score_input);
      if (kk == 0)
      {
        printf("need track the right position of kk\n");
        exit(1);
      }
      fprintf(fw, "%d\t%d\t%d\t%ld\t%ld\t%.2f\t%s\t%s\t%s\t%s\n",
        i + 1,
        kk,
        kk + 1,
        gene[kk - 1]->end,
        gene[kk]->start,
        expression_score_input[kk - 1],
        check_heg(kk, 0),
        check_nap(kk, 0),
        check_heg(kk, 1),
        check_nap(kk, 1));
    }
  }
}

/* detect valleys on operon level*/
void detectValleyO(const std::vector<double>& expression_score_input, FILE *fw)
{
  double delta;
  int kk = 0;
  int i = 0;
  for (delta = getDelta(expression_score_input) / 10; delta < getDelta(expression_score_input) / 2; delta += 3 * getDelta(expression_score_input) / 100)
  {
    fprintf(fw, "Delta-expression_score: %.2f\n", delta);
    int emi_peaks[MAX_PEAK];
    int absorp_peaks[MAX_PEAK];
    int emi_count = 0;
    int absorp_count = 0;
    int emission_first = 0;
    if (detect_peak(expression_score_input, operonNum - 1, emi_peaks, &emi_count, MAX_PEAK, absorp_peaks, &absorp_count, MAX_PEAK, delta, emission_first))
    {
      printf("There are too many peaks.\n");
      exit(1);
    }

    for (i = 0; i < absorp_count; ++i)
    {
      if (po->IS_add_nap)
        kk = which_nap(operon[absorp_peaks[i]]->gene_id[operon[absorp_peaks[i]]->gene_num - 1], operon[absorp_peaks[i] + 1]->gene_id[0], expression_score);
      if (kk == 0 || !po->IS_add_nap)
        kk = which_smallest(operon[absorp_peaks[i]]->gene_id[operon[absorp_peaks[i]]->gene_num - 1], operon[absorp_peaks[i] + 1]->gene_id[0], expression_score);
      fprintf(fw, "%d\t%d\t%d\t%ld\t%ld\t%.2f\t%s\t%s\t%s\t%s\n",
        i + 1,
        kk,
        kk + 1,
        gene[kk - 1]->end,
        gene[kk]->start,
        expression_score[kk - 1],
        check_heg(kk, 0),
        check_nap(kk, 0),
        check_heg(kk, 1),
        check_nap(kk, 1));
    }
  }
}
