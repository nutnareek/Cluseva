using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//using Cluseva;


namespace Cluseva
{
    public static class ClusterEvaluation
    {

        #region ARI
        // Implemented based on the following publication:
        // Comparing partitions
        // by Lawrence HubertPhipps Arabie
        // Accessed via https://link.springer.com/article/10.1007/BF01908075


        /// <summary>
        /// 
        /// </summary>
        /// <param name="lblSet1"></param>
        /// <param name="lblSet2"></param>
        /// <param name="numGroupSet1"></param>
        /// <param name="numGroupSet2"></param>
        /// <returns></returns>
        public static double AdjustedRandIndex(int[] lblSet1, int[] lblSet2, int numGroupSet1, int numGroupSet2)
        {
            int[][] Contingency = new int[numGroupSet2][];
            for (int i = 0; i < numGroupSet2; i++) Contingency[i] = new int[numGroupSet1];
            for (int i = 0; i < lblSet2.Length; i++)
            {
                Contingency[lblSet2[i]][lblSet1[i]] += 1;
            }

            int[] rowSums = new int[Contingency.Length];
            int[] colSums = new int[Contingency[0].Length];
            for (int i = 0; i<Contingency.Length; i++)
            {
                for(int j=0; j<Contingency[i].Length; j++)
                {
                    rowSums[i] += Contingency[i][j];
                    colSums[j] += Contingency[i][j];
                }
            }

            double sumN = 0, sumA = 0, sumB = 0;
            for (int i = 0; i < Contingency.Length; i++)
            {
                for (int j = 0; j < Contingency[i].Length; j++)
                {
                    sumN += BinomialCoefficient(Contingency[i][j], 2);
                }
            }

            for (int i = 0; i < rowSums.Length; i++)
            {
                sumA += BinomialCoefficient(rowSums[i], 2);
            }

            for (int i = 0; i < colSums.Length; i++)
            {
                sumB += BinomialCoefficient(colSums[i], 2);
            }

            return (sumN - (sumA * sumB / BinomialCoefficient(lblSet2.Length, 2))) / (((sumA + sumB) / 2.0) - (sumA * sumB / BinomialCoefficient(lblSet2.Length, 2)));
        }


        #endregion

        #region Silhouette Values
        // Implemented based on the equations in the following publication:
        // Silhouettes: A graphical aid to the interpretation and validation of cluster analysis 
        // by Peter J.Rousseeuw
        // Accessed via https://doi.org/10.1016/0377-0427(87)90125-7

        /// <summary>
        /// 
        /// </summary>
        /// <param name="clusters"></param>
        /// <param name="average"></param>
        /// <returns></returns>
        public static double[][] Silhouette(ClusterItem[][] clusters, out double average)
        {

            // calculate dissimilarity 

            double a;
            double b;
            double[][] s = new double[clusters.Length][];
            double[] interDissimilar = null;
            int[] otherClusters = null;
            int numGroup = clusters.Length;
            for (int i = 0; i < clusters.Length; i++) // i_th clustert
            {
                s[i] = new double[clusters[i].Length];
                List<int> others = new List<int>();
                for (int j = 0; j < numGroup; j++)
                    if (i != j) others.Add(j);

                // Index of other clusters's legend
                otherClusters = others.ToArray();

                for (int j = 0; j < clusters[i].Length; j++) // j_th item in cluster i_th
                {
                    if (clusters[i].Length == 1)
                    {
                        s[i][j] = 0;
                    }
                    else
                    {
                        //within cluster
                        a = AverageDissimilar(clusters[i][j].InputVector, clusters[i], new Distance());

                        // between cluster
                        interDissimilar = new double[otherClusters.Length];
                        for (int k = 0; k < otherClusters.Length; k++)
                        {
                            interDissimilar[k] = AverageInterDissimilar(clusters[i][j].InputVector, clusters[otherClusters[k]], new Distance());
                        }

                        b = interDissimilar.Min();
                        s[i][j] = SilhouetteValue(a, b);
                        clusters[i][j].SilhouetteValue = s[i][j];
                    }
                }

                // Reorder silhouette values in each cluster in descending order
                s[i] = s[i].OrderByDescending(x => x).ToArray();
            }

            average = 0;
            int totalSamples = 0;

            if (s != null)
            {
                for (int i = 0; i < s.Length; i++)
                {
                    totalSamples += s[i].Length;
                    average += s[i].Sum();
                }
                average = average / totalSamples;
            }
            else
                average = 0;
            return s;

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="labels"></param>
        /// <param name="inputs"></param>
        /// <param name="average"></param>
        /// <returns></returns>
        public static double[][] Silhouette(int[] labels, double[][] inputs,  out double average)
        {

            int[] legends = null;
            ClusterItem[][] clusters = GroupByCluster(inputs, labels, out legends);

            // calculate dissimilarity 
            double a;
            double b;
            double[][] s = new double[clusters.Length][];
            double[] interDissimilar = null;
            int[] otherClusters = null;

            for (int i = 0; i < clusters.Length; i++) // i_th clustert
            {
                s[i] = new double[clusters[i].Length];
                List<int> others = new List<int>();
                for (int j = 0; j < legends.Length; j++)
                    if (i != j) others.Add(j);

                // Index of other clusters's legend
                otherClusters = others.ToArray();
               
                for (int j = 0; j < clusters[i].Length; j++) // j_th item in cluster i_th
                {
                    if (clusters[i].Length == 1)
                    {
                        s[i][j] = 0;
                    }
                    else
                    {
                        //within cluster
                        a = AverageDissimilar(clusters[i][j].InputVector, clusters[i], new Distance());

                        // between cluster
                        interDissimilar = new double[otherClusters.Length];
                        for (int k = 0; k < otherClusters.Length; k++)
                        {
                            interDissimilar[k] = AverageInterDissimilar(clusters[i][j].InputVector, clusters[otherClusters[k]], new Distance());
                        }

                        b = interDissimilar.Min();
                        s[i][j] = SilhouetteValue(a, b);
                        clusters[i][j].SilhouetteValue = s[i][j];
                    }
                }

                // Reorder silhouette values in each cluster in descending order
                s[i] = s[i].OrderByDescending(x => x).ToArray();
            }

            average = 0;
            int totalSamples = 0;
            if (s != null)
            {
                for (int i = 0; i < s.Length; i++)
                {
                    totalSamples += s[i].Length;
                    average += s[i].Sum();
                }
                average = average / totalSamples;
            }
            else
                average = 0;

            return s;

        }

        private static double SilhouetteValue(double a, double b)
        {
            if (a < b) return 1 - (a / b);
            else if (a == b) return 0;
            else // a > b
                return (b / a) - 1;
        }

        private static double AverageDissimilar(double[] currInput, ClusterItem[] cluster, Distance distance)
        {
            double[] d = new double[cluster.Length];
            for (int i = 0; i < d.Length; i++)
            {
                d[i] = distance.Similarity(currInput, cluster[i].InputVector);

            }
            return d.Sum() / (d.Length - 1); // remove one element => itself
        }


        private static double AverageInterDissimilar(double[] currInput, ClusterItem[] cluster, Distance distance)
        {
            double[] d = new double[cluster.Length];
            for (int i = 0; i < d.Length; i++)
            {
                d[i] = distance.Similarity(currInput, cluster[i].InputVector);
            }
            return d.Average();
        }

        #endregion

        #region Information Criteria
        // Implemented based on equations from the following publication:
        // Feature-Space Clustering for fMRI Meta-Analysis
        //by Cyril Goutte, Lars Kai Hansen, Matthew G. Liptrot and Egill Rostrup
        // Accessed via https://www.ncbi.nlm.nih.gov/pubmed/11376501

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="inputs"></param>
        /// <param name="labels"></param>
        /// <param name="distance"></param>
        /// <param name="aic"></param>
        /// <param name="bic"></param>
        /// <param name="icl"></param>
        public static void InformationCriteria(double[][] centroids, double[][] inputs, int[] labels, Distance distance, out double aic, out double bic, out double icl)
        {
            int N = inputs.Length;      // No. of input vectors
            int Q = inputs[0].Length;
            // group by cluster label
            int[] legends = null;
            ClusterItem[][] clusters = GroupByCluster(inputs, labels, out legends);
            int K = clusters.Length;    // No. of clusters

            double likelihood = ComputeLikelihood(clusters, centroids, N, distance);

            aic = likelihood - ((K * Q) + 1);
            bic = likelihood - (((((K * Q) + 1) / 2.0) * Math.Log(N)));
            icl = bic;
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < N; i++) sum1 += Math.Log(i + ((K + 2.0) / 2.0));
            for (int k = 0; k < K; k++) { if (clusters[k] != null) for (int j = 0; j < clusters[k].Length; j++) sum2 += Math.Log(j + (3.0 / 2.0)); }
            icl = icl - sum1 + sum2;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="clusters"></param>
        /// <param name="numOfInputs"></param>
        /// <param name="distance"></param>
        /// <param name="aic"></param>
        /// <param name="bic"></param>
        /// <param name="icl"></param>
        public static void InformationCriteria(double[][] centroids, ClusterItem[][] clusters, int numOfInputs, Distance distance, out double aic, out double bic, out double icl)
        {
            int Q = clusters[0][0].InputVector.Length; // Number of parameters
            // group by cluster label

            int K = clusters.Length;    // No. of clusters
            double llk = ComputeLikelihood(clusters, centroids, numOfInputs, distance);
            aic = llk - ((K * Q) + 1);
            bic = llk - (((((K * Q) + 1) / 2.0) * Math.Log(numOfInputs)));
            icl = bic;
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < numOfInputs; i++) sum1 += Math.Log(i + ((K + 2.0) / 2.0));
            for (int k = 0; k < K; k++) { if (clusters[k] != null) for (int j = 0; j < clusters[k].Length; j++) sum2 += Math.Log(j + (3.0 / 2.0)); }
            icl = icl - sum1 + sum2;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="inputs"></param>
        /// <param name="labels"></param>
        /// <param name="distance"></param>
        /// <returns></returns>
        public static double AIC(double[][] centroids, double[][] inputs, int[] labels, Distance distance)
        {
            int N = inputs.Length;      // No. of input vectors
            int Q = inputs[0].Length;
            // group by cluster label
            int[] legends = null;
            ClusterItem[][] clusters = GroupByCluster(inputs, labels, out legends);
            int K = clusters.Length;    // No. of clusters

            return ComputeLikelihood(clusters, centroids, N, distance) - ((K * Q) + 1);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="clusters"></param>
        /// <param name="numOfInputs"></param>
        /// <param name="distance"></param>
        /// <returns></returns>
        public static double AIC(double[][] centroids, ClusterItem[][] clusters, int numOfInputs, Distance distance)
        {
            int Q = clusters[0][0].InputVector.Length; // Number of parameters
            // group by cluster label

            int K = clusters.Length;    // No. of clusters

            return ComputeLikelihood(clusters, centroids, numOfInputs, distance) - ((K * Q) + 1);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="inputs"></param>
        /// <param name="labels"></param>
        /// <param name="distance"></param>
        /// <returns></returns>
        public static double BIC(double[][] centroids, double[][] inputs, int[] labels, Distance distance)
        {
            int N = inputs.Length;      // No. of input vectors
            int Q = inputs[0].Length;
            // group by cluster label
            int[] legends = null;
            ClusterItem[][] clusters = GroupByCluster(inputs, labels, out legends);
            int K = clusters.Length;    // No. of clusters

            return ComputeLikelihood(clusters, centroids, N, distance) - (((((K * Q) + 1) / 2.0) * Math.Log(N)));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="clusters"></param>
        /// <param name="numOfInputs"></param>
        /// <param name="distance"></param>
        /// <returns></returns>
        public static double BIC(double[][] centroids, ClusterItem[][] clusters, int numOfInputs, Distance distance)
        {
            int Q = clusters[0][0].InputVector.Length; // Number of parameters
            // group by cluster label

            int K = clusters.Length;    // No. of clusters

            return ComputeLikelihood(clusters, centroids, numOfInputs, distance) - (((((K * Q) + 1) / 2.0) * Math.Log(numOfInputs)));
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="inputs"></param>
        /// <param name="labels"></param>
        /// <param name="distance"></param>
        /// <returns></returns>
        public static double ICL(double[][] centroids, double[][] inputs, int[] labels, Distance distance)
        {
            int N = inputs.Length;      // No. of input vectors
            int Q = inputs[0].Length;
            // group by cluster label
            int[] legends = null;
            ClusterItem[][] clusters = GroupByCluster(inputs, labels, out legends);
            int K = clusters.Length;    // No. of clusters

            double ICL = ComputeLikelihood(clusters, centroids, N, distance) - (((((K * Q) + 1) / 2.0) * Math.Log(N)));
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < N; i++) sum1 += Math.Log(i + ((K + 2.0) / 2.0));
            for (int k = 0; k < K; k++) { if (clusters[k] != null) for (int j = 0; j < clusters[k].Length; j++) sum2 += Math.Log(j + (3.0 / 2.0)); }
            return ICL - sum1 + sum2;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="centroids"></param>
        /// <param name="clusters"></param>
        /// <param name="numOfInputs"></param>
        /// <param name="distance"></param>
        /// <returns></returns>
        public static double ICL(double[][] centroids, ClusterItem[][] clusters, int numOfInputs, Distance distance)
        {
            int Q = clusters[0][0].InputVector.Length; // Number of parameters
            // group by cluster label

            int K = clusters.Length;    // No. of clusters

            double ICL = ComputeLikelihood(clusters, centroids, numOfInputs, distance) - (((((K * Q) + 1) / 2.0) * Math.Log(numOfInputs)));
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < numOfInputs; i++) sum1 += Math.Log(i + ((K + 2.0) / 2.0));
            for (int k = 0; k < K; k++) { if (clusters[k] != null) for (int j = 0; j < clusters[k].Length; j++) sum2 += Math.Log(j + (3.0 / 2.0)); }
            return ICL - sum1 + sum2;
        }

        private static double ComputeLikelihood(ClusterItem[][] clusters, double[][] centroids, int numSamples, Distance distance)
        {
            int N = numSamples;      // No. of input vectors

            // Calculate Within Cluster Variance
            double variance = WithinClusterVariance(clusters, centroids, N, distance);
            double[] diff = null;
            double likelihood = 0;
            for (int i = 0; i<clusters.Length; i++)
                for(int j = 0; j<clusters[i].Length; j++)
                {
                    diff = VectorSubtraction(clusters[i][j].InputVector, centroids[clusters[i][j].Label]);
                    likelihood += Math.Log(Math.Exp(-DotProduct(diff, diff) / (2.0 * variance)) / Math.Sqrt(2.0 * Math.PI * variance));
                }

            return likelihood;
        }

        public static double WithinClusterVariance(ClusterItem[][] clusters, double[][] centroids, int numSamples, Distance distance)
        {
            double var = 0;
            for (int i = 0; i < clusters.Length; i++)
            {
                for (int j = 0; j < clusters[i].Length; j++)
                {
                    var += distance.Similarity(clusters[i][j].InputVector, centroids[clusters[i][j].Label]);
                }
            }
            var = var / numSamples * 1.0;
            return var;
        }

        #endregion

        #region Prediction Strength
        /// Prediction Strength is implemented based on the following publication:
        /// Cluster Validation by Prediction Strength 
        /// by Robert TIBSHIRANI and Guenther WALTHER
        /// Can be accessed via http://statweb.stanford.edu/~gwalther/predictionstrength.pdf

        /// <summary>
        /// 
        /// </summary>
        /// <param name="K"></param>
        /// <param name="trainLabels"></param>
        /// <param name="clusters"></param>
        /// <returns></returns>
        public static double PredictionStrength(int K, int[] trainLabels, ClusterItem[][] clusters)
        {
            double[] CoM = null;

            //for (int i = 0; i < K; i++)
            //{
            //    if (items[i] != null)
            //        CoM[i] = CoMembership(items[i], trainLabels);
            //}
            return PredictionStrength(K, trainLabels, clusters, out CoM);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="K"></param>
        /// <param name="trainLabels"></param>
        /// <param name="clusters"></param>
        /// <param name="CoM"></param>
        /// <returns></returns>
        public static double PredictionStrength(int K, int[] trainLabels, ClusterItem[][] clusters, out double[] CoM)
        {
            CoM = new double[K];

            for (int i = 0; i < K; i++)
            {
                if (clusters[i] != null)
                    CoM[i] = CoMembership(clusters[i], trainLabels);
            }
            return CoM.Min();
        }

        private static double CoMembership(ClusterItem[] cluster, int[] trainLabels)
        {
            double sum = 0;
            int N = cluster.Length;

            for (int i = 0; i < cluster.Length; i++)
            {

                for (int j = 0; j < N; j++)
                {
                    if ((i != j) && (trainLabels[cluster[i].Index] == trainLabels[cluster[j].Index]))
                    {
                        sum += 1.0;
                    }

                }
            }
            double CoM = (N > 1) ? sum / (N * (N - 1.0)) * 1.0 : 0;

            return CoM;
        }
        #endregion

        #region Other Utility Methods
        public static ClusterItem[][] GroupByCluster(double[][] features, int[] c, out int[] clusterLegends)
        {
            // group by cluster label
            List<ClusterItem> items = new List<ClusterItem>();

            for (int i = 0; i < c.Length; i++)
            {
                items.Add(new ClusterItem(c[i], features[i], i));
            }

            var result = items.GroupBy(r => r.Label).OrderBy(r => r.Key);

            List<ClusterItem[]> clusters = new List<ClusterItem[]>();
            int[] legends = new int[result.Count()];
            int idx = 0;
            foreach (var group in result)
            {
                ClusterItem[] currCluster = group.ToArray();
                clusters.Add(currCluster);
                legends[idx] = group.Key;
                idx++;
            }
            clusterLegends = legends;
            return clusters.ToArray();
        }

        private static double DotProduct(double[] A, double[] B)
        {
            double product = 0;
            for (int i = 0; i < A.Length; i++)
                product += A[i] * B[i];
            return product;
        }

        private static double[] VectorSubtraction(double[] A, double[] B)
        {
            double[] result = new double[A.Length];
            for (int i = 0; i < A.Length; i++)
                result[i] = A[i] - B[i];
            return result;
        }

        private static double Factorial(int n)
        {
            double result = 1;
            for (int i = 2; i <= n; i++)
                result *= i;
            return result;
        }

        private static double BinomialCoefficient(int n, int k)
        {
            if (n == 0 || n == 1)
                return 0;
            else
                return Factorial(n) / (Factorial(n - k) * Factorial(k));
        }

        #endregion
    }



}
