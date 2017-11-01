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

        enum Previous
        {
            LEFT, ABOVE, DIAGONAL, ORIGIN
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

            int ALength = sequenceA.Sequence.Length;
            int BLength = sequenceB.Sequence.Length;
            //int[,] dp = new int[ALength, BLength];      //these sizes must be wrong. causes errors. Might have to do with "MaxCharactersToAlign" variable
            //Previous[,] previous = new Previous[ALength, BLength];

            score = 0;                                                
            alignment[0] = "";
            alignment[1] = "";

            
            

            result.Update(score,alignment[0],alignment[1]);                  // bundling your results into the right object type 
            return(result);
        }
    }
}
