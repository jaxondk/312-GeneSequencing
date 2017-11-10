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

            //Space: O(nm)
            int[,] dp = new int[ALength+1, BLength+1]; //the +1 is for the '-' row and column in dp/Backpointer matrices
            Backpointer[,] history = new Backpointer[ALength+1, BLength+1];
            if(banded)
                score = solveCostBanded(ALength, BLength, sequence, ref dp, ref history);
            else
                score = solveCost(ALength, BLength, sequence, ref dp, ref history);
            extractAlignments(ALength, BLength, sequence, ref history, ref alignment);
            result.Update(score,alignment[0],alignment[1]);                  // bundling your results into the right object type 
            return(result);
        }

        private int solveCostBanded(int ALength, int BLength, string[] sequence, ref int[,] dp, ref Backpointer[,] history)
        {
            int bandDist = 3;
            //Base Case: column 0 down 4 and row 0 across 5
            //O(1) for populating the base case
            dp[0, 0] = 0;
            history[0, 0] = Backpointer.ORIGIN;
            for (int i = 1; i <= bandDist; i++)
            {
                dp[i, 0] = i * 5;
                history[i, 0] = Backpointer.UP;
            }
            for (int j = 1; j <= bandDist; j++)
            {
                dp[0, j] = j * 5;
                history[0, j] = Backpointer.LEFT;
            }

            //Fill only the banded portion of the dp array. The bottom right cell is the optimal edit distance
            //O(n+m) 
            for (int i = 1; i <= ALength; i++)
            {
                int j = (i > bandDist) ? i - bandDist : 1;
                int distFromDiag = j - i;
                for (; Math.Abs(distFromDiag) <= bandDist && j <= BLength; j++, distFromDiag = j - i)
                {
                    int upCost = 5 + dp[i - 1, j];
                    int leftCost = 5 + dp[i, j - 1];
                    int diagCost = (sequence[0][i - 1] == sequence[1][j - 1]) ?
                        (-3 + dp[i - 1, j - 1]) : (1 + dp[i - 1, j - 1]); //If the sequences match at this position, -3 + diag is cost. Else, 1 + diag is cost. 
                                                                          //This is by Needleman-Wunscht study
                    int min = diagCost;
                    Backpointer prev = Backpointer.DIAGONAL;

                    if (distFromDiag != -(bandDist) && leftCost < min) //if you are -3 from the diag, the left box shouldn't be considered
                    {
                        min = leftCost;
                        prev = Backpointer.LEFT;
                    }
                    if (distFromDiag != bandDist && upCost < min) //if you are 3 from the diag, the up box shouldn't be considered
                    {
                        min = upCost;
                        prev = Backpointer.UP;
                    }
                    dp[i, j] = min;
                    history[i, j] = prev;
                }
            }
            return (history[ALength, BLength] == Backpointer.NULL) ? int.MaxValue : dp[ALength, BLength];
        }

        //Unrestricted: O(nm) | Banded: O(n+m)
        private int solveCost(int ALength, int BLength, string[] sequence, ref int[,] dp, ref Backpointer[,] history) //ref is needed to pass by reference
        {
            //Base Case: column 0 and row 0
            //O(n+m) for populating the base case
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
            //Unrestricted: O(nm) | Banded: O(n+m) since
            for (int i = 1; i <= ALength; i++)
            {
                for (int j = 1; j <= BLength; j++)
                {
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

        //O(n+m) - this is essentially taxi cab distance. You can go left at most m times and up at most n times
        //          and the work for each case is constant time
        private void extractAlignments(int ALength, int BLength, string[] sequence, ref Backpointer[,] history, ref string[] alignment)
        {
            int i = ALength;
            int j = BLength;
            Backpointer prev = history[i, j];
            if (prev == Backpointer.NULL)
            {
                alignment[0] = alignment[1] = "No Alignment Possible";
                return;
            }

            StringBuilder align0 = new StringBuilder(ALength);
            StringBuilder align1 = new StringBuilder(BLength);
            while (prev != Backpointer.ORIGIN)
            {
                switch(prev)
                {
                    case Backpointer.UP:
                        i -= 1;
                        align1.Insert(0, '-');
                        align0.Insert(0, sequence[0][i]); //since I decremented i, this is the correct character in sequence[0] to insert
                        break;
                    case Backpointer.LEFT:
                        j -= 1;
                        align0.Insert(0, '-');
                        align1.Insert(0, sequence[1][j]);
                        break;
                    case Backpointer.DIAGONAL:
                        i -= 1;
                        j -= 1;
                        align0.Insert(0, sequence[0][i]);
                        align1.Insert(0, sequence[1][j]);
                        break;
                }
                prev = history[i, j];
            }
            alignment[0] = align0.ToString();
            alignment[1] = align1.ToString();
        }
    }
}
