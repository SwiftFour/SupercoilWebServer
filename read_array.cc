/*
 * Author:Qin Ma <maqin@csbl.bmb.uga.edu>, Jan. 24, 2010
 * Usage: This is part of the supercoil package. Use, redistribute, modify without limitations.
 *
 * 1. Include two procedures for microarry input:
 * o read_continuous() would read a file with this format:
 * ----------------------------------------
 *              cond1    cond2   cond3
 *      gene1    3.14     -1.2     0.0
 *      gene2      nd      2.8     4.5
 * ----------------------------------------
 * values may possibly be any continuous value, e.g. log-ratio of
 * lumin intensity for two channels. The procedure then, for each
 * row, produces a distribution using method similar to outlier algorithm,
 * base on two tails of values (6%),
 * the middle part, is regarded as insignificant. This would discretize
 * the continuous value into classes. (If you want divide the data into
 * more levels, you can adjust the parameter r and q) See below.
 *
 * o read_discrete() would read a file with format like:
 * ----------------------------------------
 *              cond1    cond2   cond3
 *      gene1       1        1       0
 *      gene2       1       -1       0
 * ----------------------------------------
 * the symbols could be any integers (-32768~+32767) and represent distinct classes.
 * '0', however, will be ignored, and uncounted in the alter algorithms.
 * since they would represent a no-change class.
 *
 * 2. standard format of ptt file from NCBI
 * 3. standard format of gi.operon file from DOOR
 */

#include "read_array.h"
#include <algorithm>
#ifdef _WIN32
#include "ggets.h"

int getline(char **lineptr, size_t *n, FILE *stream)
{
  return fggets(lineptr, stream);
}
#endif
/************************************************************************/
/* Helper variables for tokenizer function */

static char *atom = NULL;
static char *atom1 = NULL;
static char delims[] = "\t\r\n";
static char delims1[] = "..";
static char delims2[] = " \r\n";
#define MAXC 100000
/* record the position of each discretized symbol in _symbols_ */
/* an unsigned short can hold all the values between 0 and USHRT_MAX inclusive. USHRT_MAX must be at least 65535*/
static int bb[USHRT_MAX];

/***********************************************************************/

/* emulate gnu gsl quantile function */
/*divide by the number of the data*/
static continuous quantile_from_sorted_data(const std::vector<continuous>& sorted_data, size_t n, double f)
{
  /*floor function returns the largest integral value less than or equal to x*/
  int i = static_cast<int>((n - 1)*f);
  continuous delta = (n - 1)*f - i;
  return (1 - delta)*sorted_data[i] + delta*sorted_data[i + 1];
}
/*divide by the value of the data*/
/*static continuous quantile_from_sorted_data_value(const continuous sorted_data[],size_t n, double f)
{
 return sorted_data[0]+f*(sorted_data[n-1]-sorted_data[0]);
}*/
/***********************************************************************/

static int charset_add(std::vector<discrete>& ar, discrete s)
{
  /*A signed short can hold all the values between SHRT_MIN and SHRT_MAX inclusive.SHRT_MIN is required to be -32767 or less,SHRT_MAX must be at least 32767*/
  int ps = s + SHRT_MAX;
  if (bb[ps] < 0)
  {
    bb[ps] = sigma;
    ar[sigma++] = s;
  }
  return bb[ps];
}

/***********************************************************************/

/* Matrix allocations (continuous and discrete 2d array) */

continuous** alloc2d(int rr, int cc)
{
  continuous** result;
  int i;
  result = new continuous*[rr];
  for (i = 0; i < rr; i++)
    result[i] = new continuous[cc];
  return result;
}

discrete** alloc2c(int rr, int cc)
{
  discrete** result;
  int i;
  result = new discrete*[rr];
  for (i = 0; i < rr; i++)
    result[i] = new discrete[cc];
  return result;
}

char** alloc2C(int rr, int cc)
{
  char** result;
  int i;
  result = new char*[rr];
  for (i = 0; i < rr; i++)
    result[i] = new char[cc];
  return result;
}

/*get gene number for PTT file*/
void get_gene_num(FILE* fp)
{
  size_t n = 0;
  char *line;
  /*omit the first three lines in ptt file*/
  getline(&line, &n, fp);
  getline(&line, &n, fp);
  getline(&line, &n, fp);
  while (getline(&line, &n, fp) >= 0)
  {
    geneNum++;
  }
  fseek(fp, 0, 0);
}

/*get operon number for .gi.operon file*/
void get_operon_num(FILE* fp)
{
  size_t n = 0;
  char *line;
  while (getline(&line, &n, fp) >= 0)
  {
    operonNum++;
  }
  fseek(fp, 0, 0);
}
/*get pathway number for pathway file*/
void get_pathway_num(FILE* fp)
{
  size_t n = 0;
  char *line;
  while (getline(&line, &n, fp) >= 0)
  {
    pathwayNum++;
  }
  fseek(fp, 0, 0);
}
/*get regulon number for regulon file*/
void get_regulon_num(FILE* fp)
{
  size_t n = 0;
  char *line;
  while (getline(&line, &n, fp) >= 0)
  {
    regulonNum++;
  }
  fseek(fp, 0, 0);
}
/*get NAP number for NAP file*/
void get_nap_num(FILE* fp)
{
  size_t n = 0;
  char *line;
  while (getline(&line, &n, fp) >= 0)
  {
    napNum++;
  }
  fseek(fp, 0, 0);
}
/*get palindromic*/
void get_palindromic_num(FILE* fp)
{
  size_t n = 0;
  char *line;
  while (getline(&line, &n, fp) >= 0)
  {
    palindromicNum++;
  }
  fseek(fp, 0, 0);
}
/*check the gene uniq in pathway*/
int is_appear_in_pathway(Pathway *p, int num, char *atom)
{
  int i = 0;
  for (i = 0; i < num; i++)
  {
    if (strcmp(atom, p->component[i]) == 0)
      return 1;
  }
  return 0;
}
/*check the pathway uniq in operon pathways*/
int is_appear_in_operon_pathway(int operon_id, int num, int current_id)
{
  int i = 0;
  for (i = 0; i < num; i++)
  {
    if (operon[operon_id - 1]->pathway_id[i] == current_id)
      return 1;
  }
  return 0;
}
/*check the gene uniq in regulon*/
int is_appear_in_regulon(Regulon *p, int num, char *atom)
{
  int i = 0;
  for (i = 0; i < num; i++)
  {
    if (strcmp(atom, p->component[i]) == 0)
      return 1;
  }
  return 0;
}
/*check the gene uniq in operon regulons*/
int is_appear_in_operon_regulon(int operon_id, int num, int current_id)
{
  int i = 0;
  for (i = 0; i < num; i++)
  {
    if (operon[operon_id - 1]->regulon_id[i] == current_id)
      return 1;
  }
  return 0;
}
/***********************************************************************/
/* Pre-read the datafile, retrieve gene labels and condition labels
 * as well as determine the matrix size
 */
void get_matrix_size(FILE* fp)
{
  /*size_t is the best type to use if you want to represent sizes of objects.
  * Using int to represent object sizes is likely to work on most modern systems, but it isn't guaranteed.
  */
  size_t n = 0;
  char *line;
  /*getline() reads an entire line, storing the address of the buffer containing the text into *line.
  *the buffer is null-terminated and includes the newline character, if a newline delimiter was found.
  */
  if (getline(&line, &n, fp) >= 0)
  {
    /*strtok function returns a pointer to the next token in str1, where str2 contains the delimiters that determine the token*/
    atom = strtok(line, delims);
    /*delete the first element in atom because the first element corresponding to description column*/
    atom = strtok(NULL, delims);
    while (atom != NULL)
    {
      /*if we do not set atom = strtok(NULL, delims), here is a infinite loop*/
      atom = strtok(NULL, delims);
      cols++;
    }
  }
  while (getline(&line, &n, fp) >= 0)
  {
    atom = strtok(line, delims);
    rows++;
  }
  /*fseed sets the position indicator associated with the stream to a new position defined by adding offset to a reference position specified by origin*/
  fseek(fp, 0, 0);
}

/*track the gi id in ptt file*/
int which_gi_in_ptt(long int pid)
{
  int i, id = 0;
  for (i = 0; i < geneNum; i++)
  {
    if (gene[i]->pid == pid)
    {
      id = i + 1;
      break;
    }
  }
  return id;
}

/*track the gi id in microarray file*/
int which_gi_in_microarray(long int pid)
{
  int i, j, id = 0;
  for (i = 0; i < geneNum; i++)
  {
    if (gene[i]->pid == pid)
    {
      for (j = 0; j < rows; j++)
      {
        if (strcmp(genes_n[j], gene[i]->synonym) == 0)
        {
          id = j + 1;
          break;
        }
      }
    }
  }
  return id;
}
/*track the Synonym id in ptt file*/
int which_synonym_in_ptt(char *synonym)
{
  int i, id = 0;
  for (i = 0; i < geneNum; i++)
  {
    if (strcmp(synonym, gene[i]->synonym) == 0)
    {
      id = i + 1;
      break;
    }
  }
  return id;
}
/*read the expression value*/
void read_expression(FILE* fp)
{
  size_t n = 0;
  char *line;
  int num = 0;
  while (getline(&line, &n, fp) >= 0)
  {
    atom = strtok(line, delims);
    atom = strtok(NULL, delims);
    expression_num[num++] = atof(atom);
  }
  fseek(fp, 0, 0);
}

/*read the AT-rich value*/
void read_AT_rich(FILE* fp)
{
  size_t n = 0;
  char *line;
  int num = 0;
  while (getline(&line, &n, fp) >= 0)
  {
    atom = strtok(line, delims);
    atom = strtok(NULL, delims);
    gene[num]->AT = atof(atom);
    atom = strtok(NULL, delims);
    gene[num++]->ATlocal = atof(atom);
  }
  fseek(fp, 0, 0);
}

/* read the NAP binding sites*/
void read_NAP(FILE* fp)
{
  size_t n = 0;
  char *line;
  nap.resize(napNum);
  NAP *nn;
  int num = 0;
  int id = 0;
  while (getline(&line, &n, fp) >= 0)
  {
    nn = new NAP();
    num++;
    nn->id = num;
    atom = strtok(line, delims);
    strcpy(nn->name, atom);
    atom = strtok(NULL, delims);
    nn->position = atoi(atom);
    id = searchwhichgene(nn->position, gene);
    if (id)
    {
      nn->IS_interGene = TRUE;
      nn->gene_id = id;
      nn->operon_id = gene[id - 1]->operon_id;
      if (gene[id - 1]->IS_operonEnd)
        nn->IS_interOpern = TRUE;
    }
    nap[num - 1] = nn;
  }
  fseek(fp, 0, 0);
}

/* load the highly expressed genes*/
void read_heg(FILE* fp)
{
  size_t n = 0;
  char *line;
  hegNum = 0;
  getline(&line, &n, fp);
  getline(&line, &n, fp);
  getline(&line, &n, fp);
  while (getline(&line, &n, fp) >= 0)
  {
    atom = strtok(line, delims);
    atom = strtok(NULL, delims);
    atom = strtok(NULL, delims);
    atom = strtok(NULL, delims);
    if (!which_gi_in_ptt(atoi(atom)))
      continue;
    gene[which_gi_in_ptt(atoi(atom)) - 1]->IS_heg = TRUE;
    hegNum++;
  }
  fseek(fp, 0, 0);
}

/* read the palindromic pattern*/
void read_palindromic(FILE* fp)
{
  size_t n = 0;
  char *line;
  palin.resize(palindromicNum);
  Palindromic *nn;
  int num = 0;
  int id = 0;
  while (getline(&line, &n, fp) >= 0)
  {
    nn = new Palindromic();
    num++;
    nn->id = num;
    atom = strtok(line, delims);
    nn->position = atoi(atom);
    palin[num - 1] = nn;
    id = searchwhichgene(nn->position, gene);
    if (id)
    {
      nn->IS_interGene = TRUE;
      nn->gene_id = id;
      nn->operon_id = gene[id - 1]->operon_id;
      if (gene[id - 1]->IS_operonEnd)
        nn->IS_interOpern = TRUE;
    }
  }
  fseek(fp, 0, 0);
}

/* load the gene information in ptt file*/
void read_ptt(FILE* fp)
{
  size_t n = 0;
  char *line;
  char location[30];
  gene.resize(geneNum);
  Gene *g;
  int num = 0;
  int i = 0;
  /*omit the first three lines in ptt file*/
  getline(&line, &n, fp);
  getline(&line, &n, fp);
  getline(&line, &n, fp);
  while (getline(&line, &n, fp) >= 0)
  {
    g = new Gene();
    num++;
    g->id = num;/*get the gene order*/
    atom = strtok(line, delims);
    strcpy(location, atom); /*get the gene location*/
    atom = strtok(NULL, delims);
    g->strand = atom[0]; /*get the gene strand*/
    atom = strtok(NULL, delims);
    g->length = atoi(atom);/*get the gene length*/
    atom = strtok(NULL, delims);
    g->pid = atoi(atom);/*get the gene pid*/
    atom = strtok(NULL, delims);
    strcpy(g->name, atom); /*get the gene name*/
    atom = strtok(NULL, delims);
    strcpy(g->synonym, atom); /*get the gene synonym*/
    atom1 = strtok(location, delims1);
    g->start = atoi(atom1); /*get the gene start*/
    atom1 = strtok(NULL, delims1);
    g->end = atoi(atom1); /*get the gene end*/
    gene[num - 1] = g;
  }
  for (i = 0; i < geneNum; i++)
  {
    gene[i]->pathway_num = 0;
  }
  fseek(fp, 0, 0);
}

/* load the operon information from gi.operon file*/
void read_operon(FILE* fp)
{
  operon.resize(operonNum);
  size_t n = 0;
  char *line;
  Operon *o;
  int num = 0;
  int num_1;
  while (getline(&line, &n, fp) >= 0)
  {
    o = new Operon();
    num++;
    o->id = num;
    num_1 = 0;
    atom = strtok(line, delims2);
    atom = strtok(NULL, delims2);
    if (which_gi_in_ptt(atoi(atom)) == 0)
    {
      atom = strtok(NULL, delims2);
      continue;
    }
    while (atom != NULL)
    {
      o->pid[num_1] = atoi(atom);
      o->gene_id[num_1++] = which_gi_in_ptt(atoi(atom));
      atom = strtok(NULL, delims2);
      gene[which_gi_in_ptt(o->pid[num_1 - 1]) - 1]->operon_id = num;
    }
    o->start = gene[which_gi_in_ptt(o->pid[0]) - 1]->start;
    o->end = gene[which_gi_in_ptt(o->pid[num_1 - 1]) - 1]->end;
    o->strand = gene[which_gi_in_ptt(o->pid[0]) - 1]->strand;
    o->gene_num = num_1;
    gene[which_gi_in_ptt(o->pid[num_1 - 1]) - 1]->IS_operonEnd = TRUE;
    operon[num - 1] = o;
    if (num == operonNum)
      break;
  }
  int i;
  for (i = 0; i < operonNum; i++)
  {
    operon[i]->pathway_num = 0;
  }
  fseek(fp, 0, 0);
}

/* load the pathway information from pathway file*/
void read_pathway(FILE* fp)
{
  pathway.resize(pathwayNum);
  size_t n = 0;
  char *line;
  Pathway *p;
  int num = 0;
  int num_1;
  int temp_id;
  while (getline(&line, &n, fp) >= 0)
  {
    p = new Pathway();
    p->component = alloc2C(200, 20);
    num++;
    p->id = num;
    num_1 = 0;
    atom = strtok(line, delims);
    strcpy(p->name, atom);
    atom = strtok(NULL, delims);
    atom1 = strtok(atom, delims2);
    while (atom1 != NULL)
    {
      if (is_appear_in_pathway(p, num_1, atom1) == 0)
      {
        strcpy(p->component[num_1++], atom1);
        atom1 = strtok(NULL, delims2);
        if (which_gi_in_ptt(atoi(p->component[num_1 - 1])) > 0)
        {
          temp_id = which_gi_in_ptt(atoi(p->component[num_1 - 1])) - 1;
          gene[temp_id]->pathway_id[(gene[temp_id]->pathway_num)++] = num;
          gene[temp_id]->IS_in_pathway = TRUE;
          operon[gene[temp_id]->operon_id - 1]->IS_in_pathway = TRUE;
          /*get pathway list which operon involved*/
          if (!is_appear_in_operon_pathway(gene[temp_id]->operon_id, operon[gene[temp_id]->operon_id - 1]->pathway_num, num))
            operon[gene[temp_id]->operon_id - 1]->pathway_id[(operon[gene[temp_id]->operon_id - 1]->pathway_num)++] = num;
        }
      }
      else
        atom1 = strtok(NULL, delims2);
    }
    p->capacity = num_1;
    pathway[num - 1] = p;
    if (num == pathwayNum)
      break;
  }
  fseek(fp, 0, 0);
}

/* load the regulon information from regulon file*/
void read_regulon(FILE* fp)
{
  regulon.resize(regulonNum);
  size_t n = 0;
  char *line;
  Regulon *p;
  int num = 0;
  int num_1;
  int temp_id;
  while (getline(&line, &n, fp) >= 0)
  {
    p = new Regulon();
    p->component = alloc2C(500, 20);
    num++;
    p->id = num;
    num_1 = 0;
    atom = strtok(line, delims);
    strcpy(p->name, atom);
    atom = strtok(NULL, delims);
    atom1 = strtok(atom, delims2);
    while (atom1 != NULL)
    {
      if (is_appear_in_regulon(p, num_1, atom1) == 0)
      {
        strcpy(p->component[num_1++], atom1);
        atom1 = strtok(NULL, delims2);
        if (which_gi_in_ptt(atoi(p->component[num_1 - 1])) > 0)
        {
          temp_id = which_gi_in_ptt(atoi(p->component[num_1 - 1])) - 1;
          gene[temp_id]->regulon_id[(gene[temp_id]->regulon_num)++] = num;
          gene[temp_id]->IS_in_regulon = TRUE;
          operon[gene[temp_id]->operon_id - 1]->IS_in_regulon = TRUE;
          /*get regulon list which operon involved*/
          if (!is_appear_in_operon_regulon(gene[temp_id]->operon_id, operon[gene[temp_id]->operon_id - 1]->regulon_num, num))
            operon[gene[temp_id]->operon_id - 1]->regulon_id[(operon[gene[temp_id]->operon_id - 1]->regulon_num)++] = num;
        }
      }
      else
        atom1 = strtok(NULL, delims2);
    }
    p->capacity = num_1;
    regulon[num - 1] = p;
    if (num == regulonNum)
      break;
  }
  fseek(fp, 0, 0);
}

/* Read in the labels on x and y, in microarray terms, genes(rows) and conditions(cols)*/
void read_labels(FILE* fp)
{
  int row = 0;
  int col;
  size_t n = 0;
  char *line;
  while (getline(&line, &n, fp) >= 0)
  {
    atom = strtok(line, delims);
    /*currently the first element in atom is the gene name of each row when row>=1, the 0 row corresponding to the line of condition names*/
    if (row >= 1)
    {
      strcpy(genes_n[row - 1], atom);
      if (which_synonym_in_ptt(atom) > 0)
        gene[which_synonym_in_ptt(atom) - 1]->expression_id = row;
      /*check if there exist a gene name equals to TFname by -T*/
      if (strcmp(atom, po->TFname) == 0)
      {
        TFindex = row - 1;
        printf("%d\n", TFindex);
      }
    }
    /*delete the first element in atom because the first element corresponding to description column*/
    atom = strtok(NULL, delims);
    col = 0;
    while (atom != NULL)
    {
      if (row == 0)
        strcpy(conds[col], atom);
      atom = strtok(NULL, delims);
      if (++col == cols) break;
    }
    if (++row == rows + 1) break;
  }
  fseek(fp, 0, 0);
}

/*read in the sub-gene list*/
void read_list(FILE* fp)
{
  int i = 0, j = 0;
  sub_genes_row = 0;
  char line[MAXC];
  while (fgets(line, MAXC, fp) != NULL)
  {
    atom = strtok(line, delims);
    strcpy(sub_genes[sub_genes_row], atom);
    sub_genes_row++;
  }

  /*update the sub_list*/
  sublist.resize(rows);
  for (i = 0; i < rows; i++)
    sublist[i] = FALSE;
  for (i = 0; i < sub_genes_row; i++)
    for (j = 0; j < rows; j++)
      if (strcmp(sub_genes[i], genes_n[j]) == 0)
        sublist[j] = TRUE;
}

/*read in the f3 value for each gene base on kernel density estimation*/
void read_density(FILE* fp)
{
  int k = 0;
  char line[MAXC];
  density.resize(rows);
  while (fgets(line, MAXC, fp) != NULL)
  {
    atom = strtok(line, delims);
    density[k] = atof(atom);
    k++;
  }
}

/* initialize data for discretization */
void init_dis()
{
  int row, col;
  /* store discretized values */
  symbols.resize(USHRT_MAX);
  /* memset sets the first num bytes of the block of memory pointed by ptr to the specified value
  * memset ( void * ptr, int value, size_t num )*/
  memset(bb, -1, USHRT_MAX * sizeof(*bb));
  /* always add an 'ignore' index so that symbols[0]==0*/
  charset_add(symbols, 0);
  /*initialize for arr_c*/
  arr_c = alloc2c(rows, cols);
  for (row = 0; row < rows; row++)
    for (col = 0; col < cols; col++)
      arr_c[row][col] = 0;
}

void read_discrete(FILE* fp)
{
  int row, col, i;
  init_dis();
  /* import data */
  size_t n = 0;
  char *line;
  row = 1;
  /* Skip first line with condition labels */
  getline(&line, &n, fp);
  /* read the discrete data from the second line */
  while (getline(&line, &n, fp) >= 0)
  {
    atom = strtok(line, delims);
    /*skip the first column*/
    atom = strtok(NULL, delims);
    col = 0;
    while (atom != NULL)
    {
      arr_c[row - 1][col] = charset_add(symbols, atoi(atom));
      atom = strtok(NULL, delims);
      if (++col == cols) break;
    }
    if (++row == rows + 1) break;
  }
  /* trim the leading spaceholder */
  printf("1: Discretized data contains %d classes with charset [ ", sigma);
  for (i = 0;i < sigma;i++)
    /*printf("%d ", symbols[i]); printf("]\n");*/
    printf("%d ", i); printf("]\n");
  fseek(fp, 0, 0);
}

void read_continuous(FILE* fp)
{
  int row, col;
  arr = alloc2d(rows, cols);
  /* import data */
  size_t n = 0;
  char *line;
  row = 1;
  /* ignore header line */
  getline(&line, &n, fp);
  while (getline(&line, &n, fp) >= 0)
  {
    atom = strtok(line, delims);
    /*skip the first column*/
    atom = strtok(NULL, delims);
    col = 0;
    while (atom != NULL)
    {
      /*we set all the aplha to ignore value 0*/
      /*Checks if parameter atom is either an uppercase or a lowercase alphabetic letter*/
      if (isalpha(*atom))
        arr[row - 1][col] = 0.0;
      else
        arr[row - 1][col] = atof(atom);
      atom = strtok(NULL, delims);
      if (++col == cols) break;
    }
    if (++row == rows + 1) break;
  }
  fseek(fp, 0, 0);
}

/***********************************************************************/

/* Discretize continuous values by revised outlier detection algorithm
 * see details in Algorithm Design section in paper
 */
discrete dis_value(float current, int divided, std::vector<float>& small, int cntl, std::vector<float>& big, int cntu)
{
  int i;
  float d_space = 1.0 / divided;
  for (i = 0; i < divided; i++)
  {
    if ((cntl > 0) && (current <= quantile_from_sorted_data(small, cntl, d_space * (i + 1))))
      return -i - 1;
    if ((cntu > 0) && (current >= quantile_from_sorted_data(big, cntu, 1.0 - d_space * (i + 1))))
      return i + 1;
  }
  return 0;
}

void discretize(const char* stream_nm)
{
  int row, col;
  std::vector<continuous> rowdata(cols);
  std::vector<float> big(cols), small(cols);
  int i, cntu, cntl;
  float f1, f2, f3, upper, lower;
  FILE *fw;
  fw = mustOpen(stream_nm, "w");
  init_dis();
  for (row = 0; row < rows; row++)
  {
    for (col = 0; col < cols; col++)
      rowdata[col] = arr[row][col];
    std::sort(rowdata.begin(), rowdata.end());
    f1 = quantile_from_sorted_data(rowdata, cols, 1 - po->QUANTILE);
    f2 = quantile_from_sorted_data(rowdata, cols, po->QUANTILE);
    if (po->IS_density)
    {
      f3 = density[row];
    }
    else
    {
      f3 = quantile_from_sorted_data(rowdata, cols, 0.5);
    }
    if ((f1 - f3) >= (f3 - f2))
    {
      upper = 2 * f3 - f2;
      lower = f2;
    }
    else
    {
      upper = f1;
      lower = 2 * f3 - f1;
    }
    cntu = 0; cntl = 0;
    for (i = 0; i < cols; i++)
    {
      if (rowdata[i] < lower)
      {
        small[cntl] = rowdata[i];
        cntl++;
      }
      if (rowdata[i] > upper)
      {
        big[cntu] = rowdata[i];
        cntu++;
      }
    }
    for (col = 0; col < cols; col++)
      arr_c[row][col] = charset_add(symbols, dis_value(arr[row][col], po->DIVIDED, small, cntl, big, cntu));
    fprintf(fw, "row %s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", genes_n[row], lower, upper, cntl, cntu);
  }
  progress("1: Discretization rules are written to %s", stream_nm);
  fclose(fw);
}

/* output the formatted matrix */
void write_imported(const char* stream_nm)
{
  int row, col;
  FILE *fw;
  fw = mustOpen(stream_nm, "w");
  fprintf(fw, "o");
  for (col = 0; col < cols; col++)
    fprintf(fw, "\t%s", conds[col]);
  fputc('\n', fw);
  for (row = 0; row < rows; row++)
  {
    fprintf(fw, "%s", genes_n[row]);
    for (col = 0; col < cols; col++)
      fprintf(fw, "\t%d", symbols[arr_c[row][col]]);
    fputc('\n', fw);
  }
  progress("2: Formatted data are written to %s", stream_nm);
  fclose(fw);
}

/***********************************************************************/
continuous get_KL(discrete *array, discrete *array_background, int a, int b)
{
  int i, j;
  std::vector<continuous> num, num_b;
  continuous IC = 0;
  num.resize(sigma);
  num_b.resize(sigma);
  for (i = 0; i < sigma; i++)
  {
    num[i] = 0;
    num_b[i] = 0;
  }
  for (i = 0;i < sigma;i++)
  {
    for (j = 0;j < a;j++)
      if (symbols[array[j]] == symbols[i])
        num[i]++;
    for (j = 0;j < b;j++)
      if (symbols[array_background[j]] == symbols[i])
        num_b[i]++;
  }
  for (i = 0;i < sigma;i++)
  {
    if (num[i] == 0) continue;
    if (num_b[i] == 0) continue;
    IC += (num[i] / a)*log2((num[i] * b) / (num_b[i] * a));
  }
  return IC;
}
/***********************************************************************/
/*get the freq for each pathway base on input microarray data*/
void get_pathway_freq()
{
  int i;
  int j;
  int k;
  int num;
  double num_col;
  std::vector<int> id;
  std::vector<int> id_e;
  for (i = 0; i < pathwayNum; i++)
  {
    id.resize(pathway[i]->capacity);
    id_e.resize(pathway[i]->capacity);
    for (j = 0; j < pathway[i]->capacity; j++)
    {
      id[j] = which_gi_in_ptt(atoi(pathway[i]->component[j]));
      id_e[j] = which_gi_in_microarray(atoi(pathway[i]->component[j]));
    }
    num_col = 0;
    for (k = 0; k < cols; k++)
    {
      num = 0;
      for (j = 0; j < pathway[i]->capacity; j++)
      {
        if (id[j] > 0 && id_e[j] > 0 && symbols[arr_c[id_e[j] - 1][k]] > 0)
        {
          /*printf ("\t\t%d\t%d\t%.2f\t%d\n",id[j],id_e[j],arr[id_e[j]-1][k],symbols[arr_c[id_e[j]-1][k]]);*/
          num++;
        }
      }
      /*printf ("\t%d\t%d\n",k,num);*/
      if (num >= pathway[i]->capacity * 0.6)
        num_col++;
    }
    pathway[i]->freq = num_col / cols;
  }
}

/***********************************************************************/
/*get the freq for each regulon base on input microarray data*/
void get_regulon_freq()
{
  int i;
  int j;
  int k;
  int num;
  double num_col;
  std::vector<int> id;
  std::vector<int> id_e;
  for (i = 0; i < regulonNum; i++)
  {
    id.resize(regulon[i]->capacity);
    id_e.resize(regulon[i]->capacity);
    for (j = 0; j < regulon[i]->capacity; j++)
    {
      id[j] = which_gi_in_ptt(atoi(regulon[i]->component[j]));
      id[j] = which_gi_in_microarray(atoi(regulon[i]->component[j]));
    }
    num_col = 0;
    for (k = 0; k < cols; k++)
    {
      num = 0;
      for (j = 0; j < regulon[i]->capacity; j++)
        if (id[j] > 0 && id_e[j] > 0 && symbols[arr_c[id_e[j] - 1][k]] > 0)
          num++;
      if (num >= regulon[i]->capacity * 0.6)
        num_col++;
    }
    regulon[i]->freq = num_col / cols;
  }
}

/***********************************************************************/
void print_gene_info(const char* fn)
{
  FILE *fw = mustOpen(fn, "w");
  int i = 0;
  int j = 0;
  for (i = 0; i < geneNum; i++)
  {
    if (gene[i]->IS_heg)
      fprintf(fw, "%d\t%ld\tHEG\t%ld\t%ld\t%c\t%s\t%s\t%d\t%d\t%s\t",
        gene[i]->id,
        gene[i]->pid,
        gene[i]->start,
        gene[i]->end,
        gene[i]->strand,
        gene[i]->name,
        gene[i]->synonym,
        gene[i]->operon_id,
        gene[i]->pathway_num,
        check_nap(i + 1, 0));
    else
      fprintf(fw, "%d\t%ld\t\t%ld\t%ld\t%c\t%s\t%s\t%d\t%d\t%s\t",
        gene[i]->id,
        gene[i]->pid,
        gene[i]->start,
        gene[i]->end,
        gene[i]->strand,
        gene[i]->name,
        gene[i]->synonym,
        gene[i]->operon_id,
        gene[i]->pathway_num,
        check_nap(i + 1, 0));
    for (j = 0; j < gene[i]->pathway_num; j++)
    {
      if (!j)
      {
        fprintf(fw, "%d", gene[i]->pathway_id[j]);
        continue;
      }
      fprintf(fw, " %d", gene[i]->pathway_id[j]);
    }
    fprintf(fw, "\t");
    for (j = 0; j < gene[i]->pathway_num; j++)
    {
      if (!j)
      {
        fprintf(fw, "%.2f", pathway[gene[i]->pathway_id[j] - 1]->freq);
        continue;
      }
      fprintf(fw, " %.2f", pathway[gene[i]->pathway_id[j] - 1]->freq);
    }
    if (po->IS_add_regulon)
    {
      fprintf(fw, "\t");
      fprintf(fw, "%d\t", gene[i]->regulon_num);
      for (j = 0; j < gene[i]->regulon_num; j++)
      {
        if (!j)
        {
          fprintf(fw, "%d", gene[i]->regulon_id[j]);
          continue;
        }
        fprintf(fw, " %d", gene[i]->regulon_id[j]);
      }
      fprintf(fw, "\t");
      for (j = 0; j < gene[i]->regulon_num; j++)
      {
        if (!j)
        {
          fprintf(fw, "%.2f", regulon[gene[i]->regulon_id[j] - 1]->freq);
          continue;
        }
        fprintf(fw, " %.2f", regulon[gene[i]->regulon_id[j] - 1]->freq);
      }
    }
    fprintf(fw, "\n");
  }
  fclose(fw);
  progress("3: Gene information are written to %s", fn);
}
