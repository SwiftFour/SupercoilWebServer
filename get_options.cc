/* Author:Qin Ma <maqin@csbl.bmb.uga.edu>, Jan. 22, 2010
 * Usage: This is part of bicluster package. Use, redistribution, modify without limitations
 * Process the options for the commandline tool
 */

/***************************************************************************/

#include "get_options.h"
#include <algorithm>
#include <string>

/***************************************************************************/
static const char USAGE[] =
"===================================================================\n\
[Usage]\n\
$ ./supercoil -i microarray -j ptt -o operon -p pathway[argument list]\n\
===================================================================\n\
[Input]\n\
-i : input file must be one of two tab-delimited formats\n\
  A) continuous data (default, use pre-set discretization (see -q and -r))\n\
     -------------------------------------\n\
     o        cond1    cond2    cond3\n\
     gene1      2.4      3.5     -2.4\n\
     gene2     -2.1      0.0      1.2\n\
     -------------------------------------\n\
  B) discrete data with arbitray classes (turn on -d)\n\
     use '0' for missing or insignificant data\n\
     -------------------------------------\n\
     o        cond1    cond2    cond3\n\
     gene1        1        2        2\n\
     gene2       -1        2        0\n\
     -------------------------------------\n\
-j : the standard ptt file\n\
-o : operon file from DOOR\n\
     -------------------------------------\n\
     1: 16127995 16127996 16127997 16127998\n\
     2: 16127999\n\
     3: 16128000\n\
     -------------------------------------\n\
-p : pathway file (pathway-gene mapping)\n\
     -------------------------------------\n\
     Phage packaging machinery       16129112 16129114\n\
     Phage tail fiber proteins       16129333 16129334 16129505 16128544\n\
     -------------------------------------\n\
-a : NAP binding sites information (optional)\n\
     -------------------------------------\n\
     fis     151\n\
     ihf     5219\n\
     fis     29361\n\
     -------------------------------------\n\
-e : HEG file (optional)\n\
     the standard ptt file format\n\
-A : AT-rich file (optional)\n\
     -------------------------------------\n\
     0       0.62    0.9\n\
     1       0.506172839506173       0.7\n\
     -------------------------------------\n\
-q : use quantile discretization for continuous data\n\
     default: 0.25 (see details in Method section in paper)\n\
-r : the number of ranks as which we treat the up(down)-regulated value\n\
     when discretization\n\
     default: 1\n\
-n : the proportion of pathway frequency base on current folding\n\
     default: .7\n\
-g : the proportion of pathway density base on current folding\n\
     default: .9\n\
-t : the proportion of expression level base on current folding\n\
     default: 60\n\
-c : the cutoff of spearman correlation coefficient\n\
     default: 0.6. [.4,.8)\n\
-S : the cutoff of cnt\n\
     default: 6\n\
-d : discrete data, where user should send their processed data\n\
     to different value classes, see above\n\
-m : the flag to descretize the continuous value considering kernel density estimation\n\
-M : the flag to identify the weight parameters for two kind of information: pathway and expression\n\
     default: FALSE\n\
===================================================================\n";

char* getCmdOption(char** begin, char** end, const std::string& option)
{
  char** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) return *itr;
  return NULL;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

/*argc is a count of the arguments supplied to the program and argc[] is an array of pointers to the strings which are those arguments-its type is array of pointer to char
 */
void get_options(int argc, char* argv[])
{
  bool is_valid = TRUE;

  po = new Prog_options();

  /*The getopt function gets the next option argument from the argument list specified by the argv and argc arguments.
  *Normally these values come directly from the arguments received by main
  */
  /*An option character in this string can be followed by a colon (:) to indicate that it takes a required argument.
  *If an option character is followed by two colons (::), its argument is optional
  *if an option character is followed by no colons, it does not need argument
  */

  /*optarg is set by getopt to point at the value of the option argument, for those options that accept arguments*/
  po->FN = getCmdOption(argv, argv + argc, "-i");
  po->TN = getCmdOption(argv, argv + argc, "-j");
  po->ON = getCmdOption(argv, argv + argc, "-o");
  po->PN = getCmdOption(argv, argv + argc, "-p");
  //po->BN = getCmdOption(argv, argv + argc, "-b");
  po->IS_add_nap = cmdOptionExists(argv, argv + argc, "-a"); po->AN = getCmdOption(argv, argv + argc, "-a");
  po->IS_add_heg = cmdOptionExists(argv, argv + argc, "-e"); po->EN = getCmdOption(argv, argv + argc, "-e");
  po->IS_add_regulon = cmdOptionExists(argv, argv + argc, "-R"); po->RN = getCmdOption(argv, argv + argc, "-R");
  po->IS_expression = cmdOptionExists(argv, argv + argc, "-D"); po->DN = getCmdOption(argv, argv + argc, "-D");
  po->IS_AT_rich = cmdOptionExists(argv, argv + argc, "-A"); po->CN = getCmdOption(argv, argv + argc, "-A");
  po->IS_palindromic = cmdOptionExists(argv, argv + argc, "-P"); po->GN = getCmdOption(argv, argv + argc, "-P");

  /*atof can convert string to double*/
  po->QUANTILE = cmdOptionExists(argv, argv + argc, "-q") ? std::atof(getCmdOption(argv, argv + argc, "-q")) : 0.25;
  /*atoi can convert string to integer*/
  po->DIVIDED = cmdOptionExists(argv, argv + argc, "-r") ? std::atoi(getCmdOption(argv, argv + argc, "-r")) : 1;
  po->IS_DISCRETE = cmdOptionExists(argv, argv + argc, "-d");
  po->IS_SWITCH = cmdOptionExists(argv, argv + argc, "-s");
  po->FILTER = cmdOptionExists(argv, argv + argc, "-f") ? std::atof(getCmdOption(argv, argv + argc, "-f")) : 0.25;
  po->COL_WIDTH = cmdOptionExists(argv, argv + argc, "-k") ? std::atoi(getCmdOption(argv, argv + argc, "-k")) : 5;
  po->SPEARMAN = cmdOptionExists(argv, argv + argc, "-c") ? std::atof(getCmdOption(argv, argv + argc, "-c")) : 0.6;
  /*puts writes the C string pointed by str to stdout and appends a newline character ('\n')*/
  /*exit(0) causes the program to exit with a successful termination
  *break is normally used to jump to the end of the current block of code
  *exit is normally used to shut down the current process
  */
  po->CNT = cmdOptionExists(argv, argv + argc, "-S") ? std::atoi(getCmdOption(argv, argv + argc, "-S")) : 6;
  po->IS_list = cmdOptionExists(argv, argv + argc, "-l"); po->LN = getCmdOption(argv, argv + argc, "-l");
  po->IS_density = cmdOptionExists(argv, argv + argc, "-m"); po->MN = getCmdOption(argv, argv + argc, "-m");
  po->parameter = cmdOptionExists(argv, argv + argc, "-n") ? std::atof(getCmdOption(argv, argv + argc, "-n")) : 1;
  po->parameter1 = cmdOptionExists(argv, argv + argc, "-g") ? std::atof(getCmdOption(argv, argv + argc, "-g")) : 1;
  po->parameter3 = cmdOptionExists(argv, argv + argc, "-N") ? std::atof(getCmdOption(argv, argv + argc, "-N")) : 0;
  po->parameter4 = cmdOptionExists(argv, argv + argc, "-G") ? std::atof(getCmdOption(argv, argv + argc, "-G")) : 0;
  po->parameter2 = cmdOptionExists(argv, argv + argc, "-t") ? std::atof(getCmdOption(argv, argv + argc, "-t")) : 1;
  po->uplimit = cmdOptionExists(argv, argv + argc, "-U") ? std::atoi(getCmdOption(argv, argv + argc, "-U")) : 110000;
  if (cmdOptionExists(argv, argv + argc, "-h")) { puts(USAGE); exit(0); }
  po->High = cmdOptionExists(argv, argv + argc, "-H");
  po->dynamicParameter = cmdOptionExists(argv, argv + argc, "-M");
  /* basic sanity check */
  if (is_valid && po->FN[0] == ' ')
  {
    /*errAbort("You must specify input file (-i)");*/
    puts(USAGE);
    exit(0);
  }
  if (is_valid)
    po->FF = mustOpen(po->FN, "r");
  if (is_valid)
    po->FT = mustOpen(po->TN, "r");
  if (is_valid)
    po->FO = mustOpen(po->ON, "r");
  if (is_valid)
    po->FP = mustOpen(po->PN, "r");
  if (po->IS_add_nap)
    po->FA = mustOpen(po->AN, "r");
  if (po->IS_add_heg)
    po->FE = mustOpen(po->EN, "r");
  if (po->IS_add_regulon)
    po->FR = mustOpen(po->RN, "r");
  if (po->IS_expression)
    po->FD = mustOpen(po->DN, "r");
  if (po->IS_AT_rich)
    po->FC = mustOpen(po->CN, "r");
  if (po->IS_palindromic)
    po->FG = mustOpen(po->GN, "r");
  if (po->IS_SWITCH)
  {
    po->IS_DISCRETE = TRUE;
    po->FB = mustOpen(po->BN, "r");
  }
  if (po->IS_list)
    po->FL = mustOpen(po->LN, "r");
  if (po->IS_density)
    po->FM = mustOpen(po->MN, "r");

  /* option value range check */
  if ((po->QUANTILE > .5) || (po->QUANTILE <= 0))
  {
    err("-q quantile discretization should be (0,.5]");
    is_valid = FALSE;
  }
  if (po->IS_SWITCH)
  {
    po->IS_DISCRETE = TRUE;
    po->FB = mustOpen(po->BN, "r");
  }
  if (po->IS_list)
    po->FL = mustOpen(po->LN, "r");
  if (po->IS_density)
    po->FM = mustOpen(po->MN, "r");

  /* option value range check */
  if ((po->QUANTILE > .5) || (po->QUANTILE <= 0))
  {
    err("-q quantile discretization should be (0,.5]");
    is_valid = FALSE;
  }
  if (!is_valid)
    errAbort("Type -h to view possible options");

}
/***************************************************************************/

