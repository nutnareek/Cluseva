using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Cluseva
{
    public class Distance
    {
        public enum DistanceMetric
        {
            SquaredEuclidean,
            Euclidean
        }

        public DistanceMetric Metric { get; set; }

        public Distance()
        {
            this.Metric = DistanceMetric.SquaredEuclidean;
        }

        public Distance(DistanceMetric metric)
        {
            this.Metric = metric;
        }

        public double Similarity(double[] A, double[] B)
        {
            double result = 0;
            // Default is square euclidean distance
            switch (Metric)
            {
                case DistanceMetric.Euclidean:
                    for (int i = 0; i < A.Length; i++)
                        result += Math.Pow(A[i] - B[i], 2);
                    result = Math.Sqrt(result);
                    break;
                case DistanceMetric.SquaredEuclidean:
                default:
                    for (int i = 0; i < A.Length; i++)
                        result += Math.Pow(A[i] - B[i], 2);
                    break;

            }

            return result;
        }
    }
}
