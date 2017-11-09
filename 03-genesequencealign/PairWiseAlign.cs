using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab
{
    class PairWiseAlign
    {
        int MaxCharactersToAlign;

        public PairWiseAlign()
        {
            // Default is to align only 5000 characters in each sequence.
            this.MaxCharactersToAlign = 5000;
        }

        public PairWiseAlign(int len)
        {
            // Alternatively, we can use a different length; typically used with the banded option checked.
            this.MaxCharactersToAlign = len;
        }

        enum Backpointer
        {
            NULL, LEFT, UP, DIAGONAL, ORIGIN
        };

        /**
        This is the function you implement.
        <param name="sequenceA">the first sequence</param>
        <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
        <param name="banded">true if alignment should be band limited.</param>
        <returns>the alignment score and the alignment (in a Result object) for sequenceA and sequenceB.  The calling function places the result in the dispay appropriately.
        **/
        public ResultTable.Result Align_And_Extract(GeneSequence sequenceA, GeneSequence sequenceB, bool banded)
        {
            ResultTable.Result result = new ResultTable.Result();
            int score;                                                       // place your computed alignment score here
            string[] alignment = new string[2];                             // place your two computed alignments here

            string[] sequence = new string[2];
            sequence[0] = sequenceA.Sequence;
            sequence[1] = sequenceB.Sequence;
            int ALength = (sequenceA.Sequence.Length > MaxCharactersToAlign) ? MaxCharactersToAlign : sequenceA.Sequence.Length; 
            int BLength = (sequenceB.Sequence.Length > MaxCharactersToAlign) ? MaxCharactersToAlign : sequenceB.Sequence.Length;
            int[,] dp = new int[ALength+1, BLength+1]; //the +1 is for the '-' row and column in dp/Backpointer matrices
            Backpointer[,] history = new Backpointer[ALength+1, BLength+1];

            score = solveCost(ALength, BLength, sequence, banded, ref dp, ref history);
            extractAlignments(ALength, BLength, ref history, ref alignment);
            result.Update(score,alignment[0],alignment[1]);                  // bundling your results into the right object type 
            return(result);
        }

        private int solveCost(int ALength, int BLength, string[] sequence, bool banded,
            ref int[,] dp, ref Backpointer[,] history) //ref is needed to pass by reference
        {
            //Base Case: column 0 and row 0
            dp[0, 0] = 0;
            history[0, 0] = Backpointer.ORIGIN;
            for (int i = 1; i <= ALength; i++)
            {
                dp[i, 0] = i * 5;
                history[i, 0] = Backpointer.UP;
            }
            for (int j = 1; j <= BLength; j++)
            {
                dp[0, j] = j * 5;
                history[0, j] = Backpointer.LEFT;
            }

            //Fill dp array. The bottom right cell is the optimal edit distance
            for (int i = 1; i <= ALength; i++)
            {
                for (int j = 1; j <= BLength; j++)
                {
                    int bandDist = Math.Abs(i - j);
                    if(banded && bandDist > 3)
                    {
                        dp[i, j] = int.MaxValue - 5; //5 is the most that will be added to this, and we don't want overflow
                        //history doesn't need to be set since it's initialized to NULL
                        continue;
                    }
                    int upCost = 5 + dp[i - 1, j];
                    int leftCost = 5 + dp[i, j - 1];
                    int diagCost = (sequence[0][i - 1] == sequence[1][j - 1]) ? 
                        (-3 + dp[i - 1, j - 1]) : (1 + dp[i - 1, j - 1]); //If the sequences match at this position, -3 + diag is cost. Else, 1 + diag is cost. 
                                                                            //This is by Needleman-Wunscht study
                    int min = diagCost;
                    Backpointer prev = Backpointer.DIAGONAL;

                    if(leftCost < min)
                    {
                        min = leftCost;
                        prev = Backpointer.LEFT;
                    }
                    if (upCost < min)
                    {
                        min = upCost;
                        prev = Backpointer.UP;
                    }
                    dp[i, j] = min;
                    history[i, j] = prev;
                }
            }
            return dp[ALength, BLength];
        }

        private void extractAlignments(int Alength, int BLength, ref Backpointer[,] previous, ref string[] alignment)
        {
            alignment[0] = "";
            alignment[1] = "";
        }
    }
}
