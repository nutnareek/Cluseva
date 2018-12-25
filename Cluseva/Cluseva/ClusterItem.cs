using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Cluseva
{
    public class ClusterItem
    {
        public int Index { set; get; }
        public int Label { set; get; }
        public double[] InputVector { get; set; }
        //public double Distance { get; set; }
        public double SilhouetteValue { get; set; }

        public ClusterItem() { }


        public ClusterItem(int label, double[] inputVector, int index)
        {
            Label = label;
            InputVector = inputVector;
            Index = index;
        }


    }

}
