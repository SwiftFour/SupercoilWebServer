#include "dynamic_programming.h"
#include <vector>
#include <algorithm>

/*void DyProgramming (const char* fn)*/
void DyProgramming(FILE *fw)
{
  /*FILE *fw = mustOpen(fn, "w");*/
  std::vector<Operon*> DpNode;
  std::vector<Operon*> DpNode_cp;
  Operon *temp;
  int i;
  int g;
  int c;
  int h;
  int id = 0;
  int circle;
  int start = lowlimit;
  int end = po->uplimit;
  int startoperon;
  int endoperon;
  double minedgevalue;
  int poolsum;
  solve *pollsolve[1024];
  solve *temp_solve;
  operonInPathwayNum = 0;

  /*get the qualified operon number*/
  if (po->IS_add_regulon)
  {
    for (i = 0; i < operonNum; i++)
      if (operon[i]->IS_in_pathway || operon[i]->IS_in_pathway)
        operonInPathwayNum++;
  }
  else
  {
    for (i = 0; i < operonNum; i++)
      if (operon[i]->IS_in_pathway)
        operonInPathwayNum++;
  }
  /*printf ("%d\n",operonInPathwayNum);*/

  /*get the dynamic nodes base on operon structure*/
  DpNode.resize(operonInPathwayNum);
  DpNode_cp.resize(operonInPathwayNum);
  DPsolve.resize(100);

  if (po->IS_add_regulon)
  {
    for (i = 0; i < operonNum; i++)
    {
      if (operon[i]->IS_in_pathway || operon[i]->IS_in_regulon)
      {
        temp = operon[i];
        temp->solvesum = 0;
        temp->edgevalue[0] = -2;
        DpNode[id] = temp;
        DpNode_cp[id] = temp;
        id++;
      }
    }
  }
  else
  {
    for (i = 0; i < operonNum; i++)
    {
      if (operon[i]->IS_in_pathway)
      {
        temp = operon[i];
        temp->solvesum = 0;
        temp->edgevalue[0] = -2;
        DpNode[id] = temp;
        DpNode_cp[id] = temp;
        id++;
      }
    }
  }

  /*dynamic programming in a circle genome*/
  int circleNum = 0;
  DPSOLVE *dpsolve;
  /*for(circle=0; DpNode_cp[circle]->start< po->uplimit; circle++)*/
  for (circle = 0; circle < 1; circle++)
  {
    dpsolve = new DPSOLVE();
    dpsolve->heg = alloc2C(1000, 100);
    dpsolve->nap = alloc2C(1000, 10000);
    dpsolve->palin = alloc2C(1000, 100);
    dpsolve->heg_s = alloc2C(1000, 100);
    dpsolve->nap_s = alloc2C(1000, 10000);
    dpsolve->palin_s = alloc2C(1000, 100);

    if (circle > 0)
    {
      DpNode.resize(operonInPathwayNum);
      start = lowlimit + DpNode_cp[circle]->start;
      end = po->uplimit + DpNode_cp[circle]->start;
      for (i = 0; i < operonInPathwayNum - circle; i++)
        DpNode[i] = DpNode_cp[i + circle];
      for (i = operonInPathwayNum - circle; i < operonInPathwayNum; i++)
      {
        DpNode[i] = DpNode_cp[i + circle - operonInPathwayNum];
        DpNode[i]->startnew = DpNode[i]->start + DpNode_cp[operonInPathwayNum - 1]->end;
        DpNode[i]->endnew += DpNode_cp[operonInPathwayNum - 1]->end;
      }
    }

    /*reset the solvesum and edgevalue*/
    for (i = 0; i < operonInPathwayNum; i++)
    {
      DpNode[i]->solvesum = 0;
      DpNode[i]->edgevalue[0] = -2;
      DpNode[i]->precurrentfold[0] = 1;
      if (!DpNode[i]->startnew)
        DpNode[i]->startnew = DpNode[i]->start;
      if (!DpNode[i]->endnew)
        DpNode[i]->endnew = DpNode[i]->end;
    }

    /*step I in dynamic programming*/
    startoperon = searchwhichoperon(start, DpNode);
    endoperon = searchwhichoperon(end, DpNode) - 1;

    int end_e;
    for (i = startoperon - 1; i < endoperon; i++)
    {
      DpNode[i]->solvesum++;
      /*start node*/
      DpNode[i]->precurrentfold[0] = 1;
      DpNode[i]->previoussolve[0] = 0;
      DpNode[i]->boundaryNum[0] = 2;
      int start_e = operon[DpNode[0]->id - 1]->gene_id[0];
      end_e = operon[DpNode[i]->id - 1]->gene_id[DpNode[i]->gene_num - 1];
      DpNode[i]->currentExpression[0] = (expression_num[start_e - 1] + expression_num[start_e] + expression_num[end_e - 1] + expression_num[end_e - 2] + expression_num[end_e]) / 5;

      /*calculate the ob_function_pathway*/
      if (po->IS_add_regulon)
        DpNode[i]->currentpwaysum[0] = ob_function_pathway(1, i + 1, DpNode) + ob_function_regulon(1, i + 1, DpNode);
      else
        DpNode[i]->currentpwaysum[0] = ob_function_pathway(1, i + 1, DpNode);

      /*DpNode[i]->edgevalue[0] = DpNode[i]->currentpwaysum[0];*/
      DpNode[i]->edgevalue[0] = DpNode[i]->currentpwaysum[0] + po->parameter2*DpNode[i]->currentExpression[0];
    }

    /*step II in dynamic programming*/
    start = start + lowlimit;
    for (g = searchwhichnewoperon(start, DpNode) - 1; g < operonInPathwayNum; g++)
    {
      poolsum = 0;
      int lastoperon = searchwhichnewoperon((DpNode[g]->startnew - po->uplimit) < 0 ? 0 : (DpNode[g]->startnew - po->uplimit), DpNode);
      if (lastoperon >= g)
        lastoperon = g - 1;
      /*printf ("%d\t%d\t%ld\t%ld\t%ld\t%d\t%.2f\n", g,lastoperon, DpNode[g]->start, DpNode[g]->startnew, DpNode[g]->startnew-po->uplimit, searchwhichnewoperon((DpNode[g]->startnew-po->uplimit), DpNode), DpNode[searchwhichnewoperon((DpNode[g]->startnew-po->uplimit), DpNode)]->edgevalue[0]);*/
      for (c = lastoperon; c < searchwhichnewoperon(DpNode[g]->startnew - lowlimit, DpNode) + 1;c++)
      {
        if (DpNode[c]->edgevalue[0] != -2)
        {
          temp_solve = new solve();
          temp_solve->precurrentfold = c;
          temp_solve->boundaryNum = DpNode[c]->boundaryNum[0] + 1;
          end_e = operon[DpNode[g]->id - 1]->gene_id[DpNode[g]->gene_num - 1];
          temp_solve->currentExpression = ((DpNode[c]->currentExpression[0] * DpNode[c]->boundaryNum[0]) + expression_num[end_e - 1] + expression_num[end_e - 2] + expression_num[end_e]) / (temp_solve->boundaryNum + 2);
          temp_solve->previoussolve = 1;
          if (po->IS_add_regulon)
            temp_solve->currentpwaysum = ob_function_pathway(c + 1, g + 1, DpNode) + ob_function_regulon(c + 1, g + 1, DpNode);
          else
            temp_solve->currentpwaysum = ob_function_pathway(c + 1, g + 1, DpNode);
          temp_solve->edgevalue = DpNode[c]->edgevalue[0] + temp_solve->currentpwaysum + po->parameter2*(temp_solve->currentExpression - DpNode[c]->currentExpression[0]);
          pollsolve[poolsum] = temp_solve;
          poolsum++;
        }
      }
      /*get the minimal edge value*/
      minedgevalue = pollsolve[0]->edgevalue;
      for (h = 0;h < poolsum;h++)
        if (pollsolve[h]->edgevalue <= minedgevalue)
          minedgevalue = pollsolve[h]->edgevalue;
      for (h = 0;h < poolsum;h++)
      {
        if (pollsolve[h]->edgevalue == minedgevalue)
        {
          DpNode[g]->solvesum++;
          DpNode[g]->edgevalue[DpNode[g]->solvesum - 1] = pollsolve[h]->edgevalue;
          DpNode[g]->precurrentfold[DpNode[g]->solvesum - 1] = pollsolve[h]->precurrentfold;
          DpNode[g]->boundaryNum[DpNode[g]->solvesum - 1] = pollsolve[h]->boundaryNum;
          DpNode[g]->currentExpression[DpNode[g]->solvesum - 1] = pollsolve[h]->currentExpression;
          DpNode[g]->previoussolve[DpNode[g]->solvesum - 1] = pollsolve[h]->previoussolve;
          DpNode[g]->currentpwaysum[DpNode[g]->solvesum - 1] = pollsolve[h]->currentpwaysum;
        }
      }
      /*printf ("g: %d\t%d\t%.2f\n",g, DpNode[g]->precurrentfold[0],DpNode[g]->edgevalue[0]);*/
    }
    /*output the results*/
    dpsolve->obvalue = DpNode[operonInPathwayNum - 1]->edgevalue[0];
    dpsolve->score = static_cast<int>(1000 * DpNode[operonInPathwayNum - 1]->edgevalue[0]);
    dpsolve->startOperonId = circle + 1;

    i = operonInPathwayNum - 1;
    double ii = 0, jj = 0;
    int gg, kk = 0, boundaryNum = 0;
    while (DpNode[i]->precurrentfold[0] != 1)
    {
      kk = 0;
      if (i == operonInPathwayNum - 1)
        g = i + 1 - operonInPathwayNum;
      else
        g = i + 1;
      if (operon[DpNode[g]->id - 1]->gene_id[0] < operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1])
        /*gg = operon[DpNode[g]->id-1]->gene_id[0]+geneNum;*/
        gg = geneNum;
      else
        gg = operon[DpNode[g]->id - 1]->gene_id[0];
      /*add palindromic pattern constrains*/
      if (po->IS_palindromic)
        kk = which_palin(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], gg, expression_num);
      if (kk == 0 && po->IS_add_nap)
        kk = which_nap(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], gg, expression_num);
      if (kk == 0)
        kk = which_smallest(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], gg, expression_num);
      /*save the results into dpsolve*/
      dpsolve->bgene[boundaryNum] = kk;
      dpsolve->expression[boundaryNum] = expression_num[kk - 1];
      strcpy(dpsolve->heg[boundaryNum], check_heg(kk, 0));
      strcpy(dpsolve->nap[boundaryNum], check_nap(kk, 0));
      strcpy(dpsolve->palin[boundaryNum], check_palin(kk, 0));
      strcpy(dpsolve->heg_s[boundaryNum], check_heg(kk, 1));
      strcpy(dpsolve->nap_s[boundaryNum], check_nap(kk, 1));
      strcpy(dpsolve->palin_s[boundaryNum], check_palin(kk, 1));

      if (po->IS_add_regulon)
        ii += ob_function_pathway(DpNode[i]->precurrentfold[0] - 1, i - 1, DpNode) + ob_function_regulon(DpNode[i]->precurrentfold[0] - 1, i - 1, DpNode);
      else
        ii += ob_function_pathway(DpNode[i]->precurrentfold[0] - 1, i - 1, DpNode);
      jj += expression_num[i - 1];
      i = DpNode[i]->precurrentfold[0];
      boundaryNum++;
    }

    /*printf ("\t%d\t%d\t%d\t%d\n",i,kk, DpNode[i]->precurrentfold[0],boundaryNum);*/
    if (po->IS_palindromic)
      kk = which_palin(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], operon[DpNode[i + 1]->id - 1]->gene_id[0], expression_num);
    if (kk == 0 && po->IS_add_nap)
      kk = which_nap(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], operon[DpNode[i + 1]->id - 1]->gene_id[0], expression_num);
    if (kk == 0)
      kk = which_smallest(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], operon[DpNode[i + 1]->id - 1]->gene_id[0], expression_num);
    dpsolve->bgene[boundaryNum] = kk;
    dpsolve->expression[boundaryNum] = expression_num[kk - 1];
    strcpy(dpsolve->heg[boundaryNum], check_heg(kk, 0));
    strcpy(dpsolve->palin[boundaryNum], check_palin(kk, 0));
    strcpy(dpsolve->nap[boundaryNum], check_nap(kk, 0));
    strcpy(dpsolve->heg_s[boundaryNum], check_heg(kk, 1));
    strcpy(dpsolve->palin_s[boundaryNum], check_palin(kk, 1));
    strcpy(dpsolve->nap_s[boundaryNum], check_nap(kk, 1));
    boundaryNum++;

    jj += expression_num[i - 1];
    if (po->IS_add_regulon)
      ii += ob_function_pathway(DpNode[i]->precurrentfold[0], i - 1, DpNode) +
      ob_function_regulon(DpNode[i]->precurrentfold[0], i - 1, DpNode) +
      po->parameter2*jj / boundaryNum;
    else
      ii += ob_function_pathway(DpNode[i]->precurrentfold[0], i - 1, DpNode) + po->parameter2*jj / boundaryNum;
    dpsolve->perturbation = ii - DpNode[operonInPathwayNum - 1]->edgevalue[0];
    dpsolve->boundaryNum = boundaryNum;
    DPsolve[circleNum] = dpsolve;
    /*printf ("%d\t%.2f\n", circle, DPsolve[circleNum]->obvalue);*/
    circleNum++;
  }
  print_results(fw, DPsolve, circleNum);

  /*fclose (fw);*/
  /*uglyTime("Final: Dynamic programming results are writen to %s", fn);*/
}

/*print out the dynamic programming results*/
void print_results(FILE *fw, std::vector<DPSOLVE*> DPsolve, int circleNum)
{
  sort_DP_list(DPsolve, circleNum);
  int i = 0;

  if (po->High)
  {
    fprintf(fw, "High confident boundary list:\n");
    for (i = 0; i < DPsolve[0]->boundaryNum; i++)
    {
      if ((strcmp(DPsolve[0]->heg_s[i], "HEG") == 0 && strcmp(DPsolve[0]->palin[i], "") != 0) ||
        (strcmp(DPsolve[0]->heg_s[i], "HEG") == 0 && strcmp(DPsolve[0]->nap[i], "") != 0) ||
        (strcmp(DPsolve[0]->palin[i], "") != 0 || strcmp(DPsolve[0]->nap[i], "") != 0))
      {
        {
          fprintf(fw, "%ld\t%ld\t%.2f\t%s\t%s\t%s\n",
            gene[DPsolve[0]->bgene[i] - 1]->end,
            gene[((DPsolve[0]->bgene[i] + 1) > geneNum ? (DPsolve[0]->bgene[i] + 1 - geneNum) : (DPsolve[0]->bgene[i] + 1)) - 1]->start,
            DPsolve[0]->expression[i],
            DPsolve[0]->heg_s[i],
            DPsolve[0]->palin[i],
            DPsolve[0]->nap[i]);
        }
      }
    }
    exit(1);
  }
  fprintf(fw, "Start from the %d operon:\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", DPsolve[0]->startOperonId, DPsolve[0]->obvalue, po->parameter, po->parameter1, po->parameter2, po->parameter3, po->parameter4, po->uplimit);
  for (i = 0; i < DPsolve[0]->boundaryNum; i++)
  {
    fprintf(fw, "%d\t%d\t%d\t%.2f\t%ld\t%ld\t%s\t%s\t%s\t%s\t%s\t%s\n",
      DPsolve[0]->startOperonId,
      DPsolve[0]->bgene[i],
      (DPsolve[0]->bgene[i] + 1) > geneNum ? (DPsolve[0]->bgene[i] + 1 - geneNum) : (DPsolve[0]->bgene[i] + 1),
      DPsolve[0]->expression[i],
      gene[DPsolve[0]->bgene[i] - 1]->end,
      gene[((DPsolve[0]->bgene[i] + 1) > geneNum ? (DPsolve[0]->bgene[i] + 1 - geneNum) : (DPsolve[0]->bgene[i] + 1)) - 1]->start,
      DPsolve[0]->heg[i],
      DPsolve[0]->palin[i],
      DPsolve[0]->nap[i],
      DPsolve[0]->heg_s[i],
      DPsolve[0]->palin_s[i],
      DPsolve[0]->nap_s[i]);
  }
  fprintf(fw, "Perturbation from the %d operon: %.2f\n", DPsolve[0]->startOperonId, DPsolve[0]->perturbation);

  /*generate random boundary*/
  fprintf(fw, "Random numbers start from %d operon:\n", DPsolve[0]->startOperonId);
  int *random, gg;
  random = uniq_integer(DPsolve[0]->boundaryNum, geneNum);
  for (gg = 0; gg < DPsolve[0]->boundaryNum; gg++)
  {
    fprintf(fw, "%d\t%d\t%d\t%.2f\t%ld\t%ld\t%s\t%s\t%s\t%s\t%s\t%s\n",
      DPsolve[0]->startOperonId,
      random[gg],
      (random[gg] + 1) > geneNum ? (random[gg] + 1 - geneNum) : (random[gg] + 1),
      expression_num[random[gg] - 1],
      gene[random[gg] - 1]->end,
      gene[(random[gg] + 1) > geneNum ? (random[gg] + 1 - geneNum) : (random[gg] + 1) - 1]->start,
      check_heg(random[gg], 0),
      check_palin(random[gg], 0),
      check_nap(random[gg], 0),
      check_heg(random[gg], 1),
      check_palin(random[gg], 1),
      check_heg(random[gg], 1));
  }
}

static int DP_cmpr(const void *a, const void *b)
/* compare function for qsort, increasing by score */
{
  return ((*(DPSOLVE **)a)->score - (*(DPSOLVE **)b)->score);
}

struct SortByX
{
  bool operator() (DPSOLVE* const & L, DPSOLVE* const & R) const { return L->score < R->score; }
};

static void sort_DP_list(std::vector<DPSOLVE*> el, int n)
{
  std::sort(el.begin(), el.begin() + n, SortByX());
}


/*search the oepron base on location id*/
int searchwhichoperon(int gauge, const std::vector<Operon*>& node)
{
  int h = 0;
  for (h = 0;h < operonInPathwayNum; h++)
  {
    if (gauge <= node[h]->start)
      return h + 1;
    else if (gauge < node[h]->end)
      return h + 1;
  }
  return 0;
}
/*search the oepron base on new location id*/
int searchwhichnewoperon(int gauge, const std::vector<Operon*>& node)
{
  int h = 0;
  for (h = 0;h < operonInPathwayNum; h++)
  {
    if (gauge <= node[h]->startnew)
      return h + 1;
    else if (gauge < node[h]->endnew)
      return h + 1;
  }
  return 0;
}
/*search the gene base on location id*/
int searchwhichgene(int gauge, const std::vector<Gene*>& node)
{
  int h = 0;
  for (h = 0;h < geneNum - 1; h++)
    if (gauge >= node[h]->end && gauge <= node[h + 1]->start)
      return h + 1;
  return 0;
}

/* caculate the objective functio base on regulon*/
float ob_function_regulon(int headindex, int rearindex, const std::vector<Operon*>& node)
{
  std::vector<int> allregulon(regulonNum);
  std::vector<int> allregulonfreq(regulonNum);
  std::vector<double> allregulondensity(regulonNum);
  int allregulonsum = 0;
  int g = 0;
  int h = 0;
  int u = 0;
  double OB = 0;
  bool flag = FALSE;

  for (g = 0;g < regulonNum;g++)
    allregulonfreq[g] = 0;

  for (h = headindex - 1; h < rearindex; h++)
  {
    for (u = 0;u < node[h]->regulon_num;u++)
    {
      flag = FALSE;
      for (g = 0; g < allregulonsum;g++)
      {
        if (allregulon[g] == node[h]->regulon_id[u])
        {
          flag = TRUE;
          allregulonfreq[g]++;
          break;
        }
      }
      if (!flag)
      {
        allregulonsum++;
        allregulon[allregulonsum - 1] = node[h]->regulon_id[u];
        allregulonfreq[allregulonsum - 1]++;
      }
    }
  }

  for (g = 0; g < allregulonsum; g++)
  {
    allregulondensity[g] = (rearindex - headindex + 1) / regulon[allregulon[g] - 1]->capacity;
    OB += (regulon[allregulon[g] - 1]->freq)*(po->parameter3 + po->parameter4*allregulondensity[g]);
  }
  return OB;
}

/* caculate the objective function base on pathway*/
float ob_function_pathway(int headindex, int rearindex, std::vector<Operon*> node)
{
  std::vector<int> allpathway(pathwayNum);
  std::vector<int> allpathwayfreq(pathwayNum);
  std::vector<double> allpathwaydensity(pathwayNum);
  int allpathwaysum = 0;
  int g = 0;
  int h = 0;
  int u = 0;
  double OB = 0;
  bool flag = FALSE;

  /*printf ("%d\t%d\t%d\t%d\n",headindex, rearindex, node[headindex-1]->id, node[rearindex-1]->id);*/
  for (g = 0;g < pathwayNum;g++)
    allpathwayfreq[g] = 0;

  for (h = headindex - 1; h < rearindex; h++)
  {
    /*printf ("%d\t%d\t%d\t%d\t%d\n", headindex, rearindex, h, node[h]->pathway_num, node[h]->id);*/
    for (u = 0;u < node[h]->pathway_num;u++)
    {
      flag = FALSE;
      for (g = 0; g < allpathwaysum;g++)
      {
        if (allpathway[g] == node[h]->pathway_id[u])
        {
          flag = TRUE;
          allpathwayfreq[g]++;
          break;
        }
      }
      if (!flag)
      {
        allpathwaysum++;
        allpathway[allpathwaysum - 1] = node[h]->pathway_id[u];
        allpathwayfreq[allpathwaysum - 1]++;
      }
      /*printf ("\t%d\t%d\t%d\t%d\n", u, allpathwaysum, allpathway[allpathwaysum-1], allpathwayfreq[allpathwaysum-1]);*/
    }
  }

  for (g = 0; g < allpathwaysum; g++)
  {
    /*allpathwaydensity[g] = (rearindex-headindex+1)/allpathwayfreq[g];*/
    allpathwaydensity[g] = (rearindex - headindex + 1) / pathway[allpathway[g] - 1]->capacity;
    /*consider the distance base on all operon not only in pathway*/
    /*allpathwaydensity[g] = (node[rearindex-1]->id-node[headindex-1]->id+1)/pathway[allpathway[g]-1]->capacity;*/
    OB += (pathway[allpathway[g] - 1]->freq)*(po->parameter + po->parameter1*allpathwaydensity[g]);
    /*printf ("%d\t%d\t%d\t%d\t%.2f\t%.2f\n",g,allpathwaysum,cols,rearindex-headindex+1, allpathwaydensity[g],OB);*/
  }
  return OB;
}

/*check whether related to HEG*/
const char *check_heg(int gene_order, int dis)
{
  int j = 0;
  int jj = 0;
  for (j = (gene_order - dis);j < (gene_order + dis + 1); j++)
  {
    jj = j;
    if (j < 1)
      jj = j + geneNum;
    else if (j > geneNum)
      jj = j - geneNum;
    if (gene[jj - 1]->IS_heg)
    {
      return "HEG";
    }
  }
  return "";
}

/*check whether include NAP binding sites after the gene*/
char *check_nap(int gene_order, int dis)
{
  int j = 0;
  int i = 0;
  int ii = 0;
  char *str = new char[1000];
  strcpy(str, "");
  for (j = 0; j < napNum; j++)
  {
    for (i = gene_order - dis;i < gene_order + dis + 1; i++)
    {
      ii = i;
      if (i < 1)
        ii = i + geneNum;
      else if (i > geneNum)
        ii = i - geneNum;
      if (ii == nap[j]->gene_id)
      {
        strcat(str, nap[j]->name);
        strcat(str, ",");
      }
    }
  }
  return str;
}

/*check whether include palindromic after the gene*/
const char *check_palin(int gene_order, int dis)
{
  int j = 0;
  int i = 0;
  int ii = 0;
  for (j = 0; j < palindromicNum; j++)
  {
    for (i = gene_order - dis;i < gene_order + dis + 1; i++)
    {
      ii = i;
      if (i < 1)
        ii = i + geneNum;
      else if (i > geneNum)
        ii = i - geneNum;
      if (ii == palin[j]->gene_id)
      {
        return "PALINDROMIC";
      }
    }
  }
  return "";
}

/*pick the smallest expression inter-genetic region*/
int which_smallest(int start, int end, const std::vector<double>& expression)
{
  int j = 0;
  int jj = 10000;
  int atRich = 1;
  int kk = 0;
  /*printf ("EXPRESSION: %d\t%d\t", start, end);*/
  for (j = start; j < end; j++)
  {
    if (!gene[j - 1]->IS_operonEnd)
      continue;
    if (expression[j - 1] < jj)
    {
      jj = expression[j - 1];
      if (jj > 3)
      {
        kk = j;
        continue;
      }
    }
    if (expression[j - 1] == jj || jj <= 3)
    {
      if (gene[j - 1]->IS_heg)
      {
        jj = expression[j - 1];
        kk = j;
      }
      else if (po->IS_AT_rich && (gene[(j == geneNum ? geneNum : 1, j)]->AT + gene[(j == geneNum ? geneNum : 1, j)]->ATlocal) / 2 < atRich)
      {
        jj = expression[j - 1];
        kk = j;
      }
      else if (expression[j - 1] < jj)
      {
        jj = expression[j - 1];
        kk = j;
      }
    }
  }
  if (kk == 0)
    kk = start;
  /*printf ("%d\n", kk);*/
  return kk;
}

/* pick the region with NAP */
int which_nap(int start, int end, const std::vector<double>& expression)
{
  int j = 0;
  int i = 0;
  int k = 0;
  double g = 10000;
  /*printf ("NAP: %d\t%d\n", start, end);*/
  for (j = start; j < end; j++)
  {
    if (!gene[j - 1]->IS_operonEnd)
      continue;
    /*printf ("\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",j,expression[j-1],gene[(j==geneNum?:1,j)]->AT,gene[(j==geneNum?:1,j)]->ATlocal,(gene[(j==geneNum?:1,j)]->AT+gene[(j==geneNum?:1,j)]->ATlocal)/2, expression[j-1]*(1-gene[(j==geneNum?:1,j)]->AT+gene[(j==geneNum?:1,j)]->ATlocal)/2);*/
    /*printf ("\t%d\t%.2f\n",j,expression[j-1]);*/
    for (i = 0; i < napNum; i++)
    {
      if (j == nap[i]->gene_id && expression[j - 1] < g)
      {
        k = j;
        g = expression[j - 1];
      }
      if (j == nap[i]->gene_id && expression[j - 1] == g)
      {
        if (gene[j - 1]->IS_heg)
        {
          g = expression[j - 1];
          k = j;
        }
      }
    }
  }
  /*printf ("%d\n", k);*/
  return k;
}

/* pick the region with palindromic */
int which_palin(int start, int end, const std::vector<double>& expression)
{
  int j = 0;
  int i = 0;
  int ii = 0;
  int k = 0;
  bool flag = FALSE;
  double g = 10000;
  for (j = start; j < end; j++)
  {
    /*printf ("%d\t%d\t%d\t%d\n", j, start, end, palindromicNum);*/
    if (!gene[j - 1]->IS_operonEnd)
      continue;
    for (i = 0; i < palindromicNum; i++)
    {
      /*printf ("%d\t%d\t%d\t%d\n", start, end, i, palin[i]->gene_id);*/
      if (j == palin[i]->gene_id && expression[j - 1] < g)
      {
        k = j;
        g = expression[j - 1];
      }
      if (j == palin[i]->gene_id && expression[j - 1] == g)
      {
        for (ii = 0; ii < napNum; ii++)
        {
          if (j == nap[ii]->gene_id)
          {
            k = j;
            flag = TRUE;
          }
        }
        if (!flag && gene[j - 1]->IS_heg)
        {
          g = expression[j - 1];
          k = j;
        }
      }
    }
  }
  return k;
}

/* generate uniq random integers*/
int *uniq_integer(int num, int NUM)
{
  std::vector<int> rand_array;

  int i = 0, rnd;
  int* vektor = new int[num];
  rand_array.resize(NUM);
  for (i = 0;i < NUM; i++)
    rand_array[i] = 0;
#ifndef _DEBUG
  srand(time(NULL));
#endif
  i = 0;
  while (i < num)
  {
    rnd = rand() % NUM + 1;
    if (rand_array[rnd - 1] == 0)
    {
      vektor[i++] = rnd;
      rand_array[rnd - 1] = 1;
    }
  }
  return vektor;
}

void DyProgramming_single(const char* fn)
{
  FILE *fw = mustOpen(fn, "w");
  std::vector<Operon*> DpNode;
  std::vector<Operon*> DpNode_cp;
  Operon *temp;
  int i = 0;
  int g = 0;
  int c = 0;
  int h = 0;
  int id = 0;
  int circle = 0;
  int start = lowlimit;
  int end = po->uplimit;
  int startoperon = 0;
  int endoperon = 0;
  double minedgevalue = 0;
  int poolsum;
  solve *pollsolve[1024];
  solve *temp_solve;
  operonInPathwayNum = 0;

  /*get the qualified operon number*/
  if (po->IS_add_regulon)
  {
    for (i = 0; i < operonNum; i++)
      if (operon[i]->IS_in_pathway || operon[i]->IS_in_pathway)
        operonInPathwayNum++;
  }
  else
  {
    for (i = 0; i < operonNum; i++)
      if (operon[i]->IS_in_pathway)
        operonInPathwayNum++;
  }
  /*printf ("%d\n",operonInPathwayNum);*/

  /*get the dynamic nodes base on operon structure*/
  DpNode.resize(operonInPathwayNum);
  DpNode_cp.resize(operonInPathwayNum);
  DPsolve.resize(100);

  if (po->IS_add_regulon)
  {
    for (i = 0; i < operonNum; i++)
    {
      if (operon[i]->IS_in_pathway || operon[i]->IS_in_regulon)
      {
        temp = operon[i];
        temp->solvesum = 0;
        temp->edgevalue[0] = -2;
        DpNode[id] = temp;
        DpNode_cp[id] = temp;
        id++;
      }
    }
  }
  else
  {
    for (i = 0; i < operonNum; i++)
    {
      if (operon[i]->IS_in_pathway)
      {
        temp = operon[i];
        temp->solvesum = 0;
        temp->edgevalue[0] = -2;
        DpNode[id] = temp;
        DpNode_cp[id] = temp;
        id++;
      }
    }
  }

  /*dynamic programming in a circle genome*/
  int circleNum = 0;
  DPSOLVE *dpsolve;
  /*for(circle=0; DpNode_cp[circle]->start< po->uplimit; circle++)*/
  for (circle = 0; circle < 1; circle++)
  {
    dpsolve = new DPSOLVE();
    dpsolve->heg = alloc2C(1000, 100);
    dpsolve->nap = alloc2C(1000, 10000);
    dpsolve->palin = alloc2C(1000, 100);
    dpsolve->heg_s = alloc2C(1000, 100);
    dpsolve->nap_s = alloc2C(1000, 10000);
    dpsolve->palin_s = alloc2C(1000, 100);

    if (circle > 0)
    {
      DpNode.resize(operonInPathwayNum);
      start = lowlimit + DpNode_cp[circle]->start;
      end = po->uplimit + DpNode_cp[circle]->start;
      for (i = 0; i < operonInPathwayNum - circle; i++)
        DpNode[i] = DpNode_cp[i + circle];
      for (i = operonInPathwayNum - circle; i < operonInPathwayNum; i++)
      {
        DpNode[i] = DpNode_cp[i + circle - operonInPathwayNum];
        DpNode[i]->startnew = DpNode[i]->start + DpNode_cp[operonInPathwayNum - 1]->end;
        DpNode[i]->endnew += DpNode_cp[operonInPathwayNum - 1]->end;
      }
    }

    /*reset the solvesum and edgevalue*/
    for (i = 0; i < operonInPathwayNum; i++)
    {
      DpNode[i]->solvesum = 0;
      DpNode[i]->edgevalue[0] = -2;
      DpNode[i]->precurrentfold[0] = 1;
      if (!DpNode[i]->startnew)
        DpNode[i]->startnew = DpNode[i]->start;
      if (!DpNode[i]->endnew)
        DpNode[i]->endnew = DpNode[i]->end;
    }

    /*step I in dynamic programming*/
    startoperon = searchwhichoperon(start, DpNode);
    endoperon = searchwhichoperon(end, DpNode) - 1;

    int start_e = 0, end_e = 0;
    for (i = startoperon - 1; i < endoperon; i++)
    {
      DpNode[i]->solvesum++;
      /*start node*/
      DpNode[i]->precurrentfold[0] = 1;
      DpNode[i]->previoussolve[0] = 0;
      DpNode[i]->boundaryNum[0] = 2;
      start_e = operon[DpNode[0]->id - 1]->gene_id[0];
      end_e = operon[DpNode[i]->id - 1]->gene_id[DpNode[i]->gene_num - 1];
      DpNode[i]->currentExpression[0] = (expression_num[start_e - 1] + expression_num[start_e] + expression_num[end_e - 1] + expression_num[end_e - 2] + expression_num[end_e]) / 5;

      /*calculate the ob_function_pathway*/
      if (po->IS_add_regulon)
        DpNode[i]->currentpwaysum[0] = ob_function_pathway(1, i + 1, DpNode) + ob_function_regulon(1, i + 1, DpNode);
      else
        DpNode[i]->currentpwaysum[0] = ob_function_pathway(1, i + 1, DpNode);

      /*DpNode[i]->edgevalue[0] = DpNode[i]->currentpwaysum[0];*/
      DpNode[i]->edgevalue[0] = DpNode[i]->currentpwaysum[0] + po->parameter2*DpNode[i]->currentExpression[0];
    }

    /*step II in dynamic programming*/
    start = start + lowlimit;
    int lastoperon = 0;
    for (g = searchwhichnewoperon(start, DpNode) - 1; g < operonInPathwayNum; g++)
    {
      poolsum = 0;
      lastoperon = searchwhichnewoperon((DpNode[g]->startnew - po->uplimit) < 0 ? 0 : (DpNode[g]->startnew - po->uplimit), DpNode);
      if (lastoperon >= g)
        lastoperon = g - 1;
      /*printf ("%d\t%d\t%ld\t%ld\t%ld\t%d\t%.2f\n", g,lastoperon, DpNode[g]->start, DpNode[g]->startnew, DpNode[g]->startnew-po->uplimit, searchwhichnewoperon((DpNode[g]->startnew-po->uplimit), DpNode), DpNode[searchwhichnewoperon((DpNode[g]->startnew-po->uplimit), DpNode)]->edgevalue[0]);*/
      for (c = lastoperon; c < searchwhichnewoperon(DpNode[g]->startnew - lowlimit, DpNode) + 1;c++)
      {
        if (DpNode[c]->edgevalue[0] != -2)
        {
          temp_solve = new solve();
          temp_solve->precurrentfold = c;
          temp_solve->boundaryNum = DpNode[c]->boundaryNum[0] + 1;
          end_e = operon[DpNode[g]->id - 1]->gene_id[DpNode[g]->gene_num - 1];
          temp_solve->currentExpression = ((DpNode[c]->currentExpression[0] * DpNode[c]->boundaryNum[0]) + expression_num[end_e - 1] + expression_num[end_e - 2] + expression_num[end_e]) / (temp_solve->boundaryNum + 2);
          temp_solve->previoussolve = 1;
          if (po->IS_add_regulon)
            temp_solve->currentpwaysum = ob_function_pathway(c + 1, g + 1, DpNode) + ob_function_regulon(c + 1, g + 1, DpNode);
          else
            temp_solve->currentpwaysum = ob_function_pathway(c + 1, g + 1, DpNode);
          temp_solve->edgevalue = DpNode[c]->edgevalue[0] + temp_solve->currentpwaysum + po->parameter2*(temp_solve->currentExpression - DpNode[c]->currentExpression[0]);
          pollsolve[poolsum] = temp_solve;
          poolsum++;
        }
      }
      /*get the minimal edge value*/
      minedgevalue = pollsolve[0]->edgevalue;
      for (h = 0;h < poolsum;h++)
        if (pollsolve[h]->edgevalue <= minedgevalue)
          minedgevalue = pollsolve[h]->edgevalue;
      for (h = 0;h < poolsum;h++)
      {
        if (pollsolve[h]->edgevalue == minedgevalue)
        {
          DpNode[g]->solvesum++;
          DpNode[g]->edgevalue[DpNode[g]->solvesum - 1] = pollsolve[h]->edgevalue;
          DpNode[g]->precurrentfold[DpNode[g]->solvesum - 1] = pollsolve[h]->precurrentfold;
          DpNode[g]->boundaryNum[DpNode[g]->solvesum - 1] = pollsolve[h]->boundaryNum;
          DpNode[g]->currentExpression[DpNode[g]->solvesum - 1] = pollsolve[h]->currentExpression;
          DpNode[g]->previoussolve[DpNode[g]->solvesum - 1] = pollsolve[h]->previoussolve;
          DpNode[g]->currentpwaysum[DpNode[g]->solvesum - 1] = pollsolve[h]->currentpwaysum;
        }
      }
      /*printf ("g: %d\t%d\t%.2f\n",g, DpNode[g]->precurrentfold[0],DpNode[g]->edgevalue[0]);*/
    }
    /*output the results*/
    dpsolve->obvalue = DpNode[operonInPathwayNum - 1]->edgevalue[0];
    dpsolve->score = static_cast<int>(1000 * DpNode[operonInPathwayNum - 1]->edgevalue[0]);
    dpsolve->startOperonId = circle + 1;

    i = operonInPathwayNum - 1;
    double ii = 0, jj = 0;
    int gg = 0, kk = 0, boundaryNum = 0;
    while (DpNode[i]->precurrentfold[0] != 1)
    {
      kk = 0;
      if (i == operonInPathwayNum - 1)
        g = i + 1 - operonInPathwayNum;
      else
        g = i + 1;
      if (operon[DpNode[g]->id - 1]->gene_id[0] < operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1])
        /*gg = operon[DpNode[g]->id-1]->gene_id[0]+geneNum;*/
        gg = geneNum;
      else
        gg = operon[DpNode[g]->id - 1]->gene_id[0];
      /*add palindromic pattern constrains*/
      if (po->IS_palindromic)
        kk = which_palin(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], gg, expression_num);
      if (kk == 0 && po->IS_add_nap)
        kk = which_nap(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], gg, expression_num);
      if (kk == 0)
        kk = which_smallest(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], gg, expression_num);
      /*save the results into dpsolve*/
      dpsolve->bgene[boundaryNum] = kk;
      dpsolve->expression[boundaryNum] = expression_num[kk - 1];
      strcpy(dpsolve->heg[boundaryNum], check_heg(kk, 0));
      strcpy(dpsolve->nap[boundaryNum], check_nap(kk, 0));
      strcpy(dpsolve->palin[boundaryNum], check_palin(kk, 0));
      strcpy(dpsolve->heg_s[boundaryNum], check_heg(kk, 1));
      strcpy(dpsolve->nap_s[boundaryNum], check_nap(kk, 1));
      strcpy(dpsolve->palin_s[boundaryNum], check_palin(kk, 1));

      if (po->IS_add_regulon)
        ii += ob_function_pathway(DpNode[i]->precurrentfold[0] - 1, i - 1, DpNode) + ob_function_regulon(DpNode[i]->precurrentfold[0] - 1, i - 1, DpNode);
      else
        ii += ob_function_pathway(DpNode[i]->precurrentfold[0] - 1, i - 1, DpNode);
      jj += expression_num[i - 1];
      i = DpNode[i]->precurrentfold[0];
      boundaryNum++;
    }

    /*printf ("\t%d\t%d\t%d\t%d\n",i,kk, DpNode[i]->precurrentfold[0],boundaryNum);*/
    if (po->IS_palindromic)
      kk = which_palin(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], operon[DpNode[i + 1]->id - 1]->gene_id[0], expression_num);
    if (kk == 0 && po->IS_add_nap)
      kk = which_nap(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], operon[DpNode[i + 1]->id - 1]->gene_id[0], expression_num);
    if (kk == 0)
      kk = which_smallest(operon[DpNode[i]->id - 1]->gene_id[operon[DpNode[i]->id - 1]->gene_num - 1], operon[DpNode[i + 1]->id - 1]->gene_id[0], expression_num);
    dpsolve->bgene[boundaryNum] = kk;
    dpsolve->expression[boundaryNum] = expression_num[kk - 1];
    strcpy(dpsolve->heg[boundaryNum], check_heg(kk, 0));
    strcpy(dpsolve->palin[boundaryNum], check_palin(kk, 0));
    strcpy(dpsolve->nap[boundaryNum], check_nap(kk, 0));
    strcpy(dpsolve->heg_s[boundaryNum], check_heg(kk, 1));
    strcpy(dpsolve->palin_s[boundaryNum], check_palin(kk, 1));
    strcpy(dpsolve->nap_s[boundaryNum], check_nap(kk, 1));
    boundaryNum++;

    jj += expression_num[i - 1];
    if (po->IS_add_regulon)
      ii += ob_function_pathway(DpNode[i]->precurrentfold[0], i - 1, DpNode) +
      ob_function_regulon(DpNode[i]->precurrentfold[0], i - 1, DpNode) +
      po->parameter2*jj / boundaryNum;
    else
      ii += ob_function_pathway(DpNode[i]->precurrentfold[0], i - 1, DpNode) + po->parameter2*jj / boundaryNum;
    dpsolve->perturbation = ii - DpNode[operonInPathwayNum - 1]->edgevalue[0];
    dpsolve->boundaryNum = boundaryNum;
    DPsolve[circleNum] = dpsolve;
    /*printf ("%d\t%.2f\n", circle, DPsolve[circleNum]->obvalue);*/
    circleNum++;
  }
  print_results(fw, DPsolve, circleNum);

  fclose(fw);
  uglyTime("Final: Dynamic programming results are writen to %s", fn);
}

