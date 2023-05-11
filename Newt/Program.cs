using MathNet.Numerics.LinearAlgebra;

namespace Newt
{
    public class Program
    {
        #region 底层函数
        static int scope = 10;

        static List<double> fig = new List<double>();

        static int iter = 100000;

        static double tol = 0.001;

        static double GetFuncNum(int scope, double[] x)
        {
            double f = 0;
            for (int i = 0; i < scope; i++)
            {
                f += Math.Pow(-x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2], 2) + Math.Pow(x[i] * x[i] - x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2], 2) + Math.Pow(x[i] * x[i] + x[i + 1] * x[i + 1] - x[i + 2] * x[i + 2], 2);
            }
            return f;
        }

        static double[] GetFuncDen(int scope, double[] x)
        {
            double[] f = new double[scope + 2];
            for (int i = 0; i < scope; i++)
            {
                f[i] += 4 * x[i] * (3 * x[i] * x[i] - x[i + 1] * x[i + 1] - x[i + 2] * x[i + 2]);
                f[i + 1] += 4 * x[i + 1] * (3 * x[i + 1] * x[i + 1] - x[i] * x[i] - x[i + 2] * x[i + 2]);
                f[i + 2] += 4 * x[i + 2] * (3 * x[i + 2] * x[i + 2] - x[i + 1] * x[i + 1] - x[i] * x[i]);
            }
            return f;
        }

        static double[,] GetFuncDen2(int scope, double[] x)
        {
            double[,] f = new double[scope + 2, scope + 2];
            for (int i = 0; i < scope; i++)
            {
                f[i, i] += 36 * x[i] * x[i] - 4 * x[i + 1] * x[i + 1] - 4 * x[i + 2] * x[i + 2];
                f[i, i + 1] += -8 * x[i] * x[i + 1];
                f[i, i + 2] += -8 * x[i] * x[i + 2];
                f[i + 1, i] += -8 * x[i] * x[i + 1];
                f[i + 1, i + 1] += 36 * x[i + 1] * x[i + 1] - 4 * x[i] * x[i] - 4 * x[i + 2] * x[i + 2];
                f[i + 1, i + 2] += -8 * x[i + 1] * x[i + 2];
                f[i + 2, i] += -8 * x[i] * x[i + 2]; ;
                f[i + 2, i + 1] += -8 * x[i + 1] * x[i + 2];
                f[i + 2, i + 2] += 36 * x[i + 2] * x[i + 2] - 4 * x[i + 1] * x[i + 1] - 4 * x[i] * x[i];
            }
            return f;
        }

        static double[] GetX(int scope)
        {
            double[] x = new double[scope + 2];
            if (scope % 3 == 1)
            {
                int count = 0;
                for (int i = 0; i < scope + 2; i++)
                {
                    switch (count)
                    {
                        case 0:
                            x[i] = 1; break;
                        case 1:
                            x[i] = 2; break;
                        default:
                            x[i] = 1; break;
                    }
                    ++count;
                    if (count == 3) count = 0;
                }
                return x;
            }
            else
            {
                Console.WriteLine("非法变量数");
                return x;
            }
        }

        static bool GetTol(int scope, double[] x)
        {
            List<double> f = new List<double>();
            double tol_cur = 0;
            double[] tol_den = GetFuncDen(scope, x);
            for (int i = 0; i < scope + 2; i++)
            {
                tol_cur += Math.Pow(tol_den[i], 2);
            }
            tol_cur = Math.Sqrt(tol_cur);
            fig.Add(tol_cur);
            return tol_cur < tol;
        }
        #endregion


        public static void Main(string[] args)
        {
            int[] count = new int[5];

            double[] x_ini1 = GetX(scope);
            double[] x11 = new double[scope + 2];
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    x11 = CalNewt(scope, x_ini1);
                }
                else
                {
                    x11 = CalNewt(scope, x11);
                }
                if (GetTol(scope, x11))
                {
                    count[0] = l; break;
                }
            }

            double[] x_ini2 = GetX(scope);
            double[] x22 = new double[scope + 2];
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    x22 = CalZNewt(scope, x_ini2);
                }
                else
                {
                    x22 = CalZNewt(scope, x22);
                }
                if (GetTol(scope, x22))
                {
                    count[1] = l; break;
                }
            }

            double[] x_ini3 = GetX(scope);
            double[] x33 = new double[scope + 2];
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    x33 = CalDFP(scope, x_ini3, l);
                }
                else
                {
                    x33 = CalDFP(scope, x33, l);
                }
                if (GetTol(scope, x33))
                {
                    count[2] = l; break;
                }
            }

            double[] x_ini4 = GetX(scope);
            double[] x44 = new double[scope + 2];
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    x44 = CalTR(scope, x_ini4, l);
                }
                else
                {
                    x44 = CalTR(scope, x44, l);
                }
                if (GetTol(scope, x44))
                {
                    count[3] = l; break;
                }
            }

            double[] x_ini5 = GetX(scope);
            double[] x55 = new double[scope + 2];
            for (int l = 0; l < iter; l++)
            {
                if (l == 0)
                {
                    x55 = CalConjGrads(scope, x_ini5, l);
                }
                else
                {
                    x55 = CalConjGrads(scope, x55, l);
                }
                if (GetTol(scope, x55))
                {
                    count[4] = l; break;
                }
            }

            foreach (var item in fig)
            {
                Console.WriteLine(item);
            }
        }

        #region newton法
        static double[] CalNewt(int scope, double[] x)
        {
            double[] den = GetFuncDen(scope, x);
            double[,] den2 = GetFuncDen2(scope, x);
            Matrix<double> hessian = Matrix<double>.Build.DenseOfArray(den2);
            Matrix<double> hessian_inverse = hessian.Inverse();
            double[,] den2_inverse = hessian_inverse.ToArray();
            double[] xm = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    xm[i] += den2_inverse[i, j] * den[j];
                }
                xm[i] = x[i] - xm[i];
            }
            return xm;
        }
        #endregion


        #region 阻尼牛顿法
        static int iterZ = 100;
        //移动倍率
        static int scopeZ = 2;
        //黄金分割率
        static double ratioZ = 0.618;
        //精度
        static double tolZ = 0.0001;
        static double[] CalZNewt(int scope, double[] x)
        {
            double[] den = GetFuncDen(scope, x);
            double[,] den2 = GetFuncDen2(scope, x);
            Matrix<double> hessian = Matrix<double>.Build.DenseOfArray(den2);
            Matrix<double> hessian_inverse = hessian.Inverse();
            double[,] den2_inverse = hessian_inverse.ToArray();
            double[] xm = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    xm[i] += den2_inverse[i, j] * den[j];
                }
            }
            return Z(xm, x);
        }

        //黄金分割法求区间
        static double[] Z(double[] xm, double[] x)
        {
            double[] a = new double[xm.Length];
            double[] b = new double[xm.Length];
            double[] x1 = new double[xm.Length];
            double[] x2 = new double[xm.Length];
            for (int i = 0; i < xm.Length; i++)
            {
                a[i] = x[i] - scopeZ * xm[i];
                b[i] = x[i];
                x1[i] = b[i] - ratioZ * (b[i] - a[i]);
                x2[i] = a[i] + ratioZ * (b[i] - a[i]);
            }
            double f1 = GetFuncNum(scope, x1);
            double f2 = GetFuncNum(scope, x2);

            for (int i = 0; i < iterZ; i++)
            {
                if (f1 > f2)
                {
                    for (int k = 0; k < xm.Length; k++)
                    {
                        a[k] = x1[k]; x1[k] = x2[k];
                        x2[k] = a[k] + ratioZ * (b[k] - a[k]);
                    }
                    f1 = f2;
                    f2 = GetFuncNum(scope, x2);
                }
                else
                {
                    for (int k = 0; k < xm.Length; k++)
                    {
                        b[k] = x2[k]; x2[k] = x1[k];
                        x1[k] = b[k] - ratioZ * (b[k] - a[k]);
                    }
                    f2 = f1;
                    f1 = GetFuncNum(scope, x1);
                }
                double temp = 0;
                for (int j = 0; j < xm.Length; j++)
                {
                    if (Math.Abs(b[j] - a[j]) >= temp)
                    {
                        temp = Math.Abs(b[j] - a[j]);
                    }
                }
                if (temp < tolZ)
                {
                    break;
                }
            }

            double[] lambda = new double[xm.Length];
            for (int i = 0; i < xm.Length; i++)
            {
                lambda[i] = 0.5 * (a[i] + b[i]);
            }
            return lambda;
        }

        #endregion

        #region DFP
        //记忆校正矩阵所需的前一次迭代的数据
        static double[] x1 = new double[scope + 2];
        static double[] den1 = new double[scope + 2];
        static double[,] aMatrix = new double[scope + 2, scope + 2];
        static double[] CalDFP(int scope, double[] x, int l)
        {
            if (l == 0)
            {

                den1 = GetFuncDen(scope, x);
                double[] xm = Z(den1, x);
                for (int i = 0; i < scope + 2; i++)
                {
                    x1[i] = x[i];
                    aMatrix[i, i] = 1;
                }
                return xm;
            }
            else
            {
                double[] den = GetFuncDen(scope, x);
                double[,] aMatrix1 = GetRevise(scope, x, den);
                double[] xm = new double[scope + 2];
                for (int i = 0; i < scope + 2; i++)
                {
                    for (int j = 0; j < scope + 2; j++)
                    {
                        xm[i] += aMatrix1[i, j] * den[j];
                    }
                }
                //记忆数据
                for (int i = 0; i < scope + 2; i++)
                {
                    x1[i] = x[i];
                    den1[i] = den[i];
                    for (int j = 0; j < scope + 2; j++)
                    {
                        aMatrix[i, j] = aMatrix1[i, j];
                    }
                }
                return Z(xm, x);
            }
        }

        static double[,] GetRevise(int scope, double[] x, double[] den)
        {
            double[,] eMatrix = new double[scope + 2, scope + 2];
            double xg = 0; double gag = 0;
            double[] gag1 = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                xg += (x[i] - x1[i]) * (den[i] - den1[i]);
                double gag2 = 0;
                for (int j = 0; j < scope + 2; j++)
                {
                    gag2 += (den[j] - den1[j]) * aMatrix[j, i];
                }
                gag1[i] = gag2;
            }
            for (int i = 0; i < scope + 2; i++)
            {
                gag += gag1[i] * (den[i] - den1[i]);
            }

            double[] ag = new double[scope + 2];
            double[] ga = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                double ag1 = 0; double ga1 = 0;
                for (int j = 0; j < scope + 2; j++)
                {
                    ag1 += aMatrix[i, j] * (den[j] - den1[j]);
                    ga1 += (den[j] - den1[j]) * aMatrix[i, j];
                }
                ag[i] = ag1;
                ga[i] = ga1;
            }

            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    eMatrix[i, j] = (x[i] - x1[i]) * (x[j] - x1[j]) / xg - ag[i] * ga[j] / gag;
                }
            }

            double[,] aMatrix1 = new double[scope + 2, scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    aMatrix1[i, j] = eMatrix[i, j] + aMatrix[i, j];
                }
            }

            return aMatrix1;
        }
        #endregion

        #region 信赖域算法

        //选取牛顿方向作为子问题的最优解方向，因为本例中B和Hessian等价

        static double delta_Max = 1.5;
        static double delta_yeta = 0.2;
        static double rou1 = 0.25;
        static double rou2 = 0.75;
        static double gamma1 = 0.25;
        static double gamma2 = 2;
        static double delta;
        static double[] CalTR(int scope, double[] x, int l)
        {
            double[] den = GetFuncDen(scope, x);
            double[,] den2 = GetFuncDen2(scope, x);

            if (l == 0)
            {
                delta = 1;
            }
            double[] xm = CalCauchyPoint(delta, den, den2);
            double xm_Norm = 0;
            for (int i = 0; i < scope + 2; i++)
            {
                xm_Norm += den[i] * den[i];
            }
            xm_Norm = Math.Sqrt(xm_Norm);

            double md = 0;
            double[] md1 = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    md1[i] += 0.5 * xm[j] * den2[j, i];
                }
                md += (md1[i] * xm[i] + den[i] * xm[i]);
            }
            md += GetFuncNum(scope, x);

            double[] x_Iter = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                x_Iter[i] = x[i] + xm[i];
            }

            double rou = (GetFuncNum(scope, x) - GetFuncNum(scope, x_Iter)) / (GetFuncNum(scope, x) - md);

            //更新信赖域
            if (rou < rou1)
            {
                delta = gamma1 * delta;
            }
            else if (rou > rou2 && xm_Norm == delta)
            {
                delta = Math.Min(gamma2 * delta, delta_Max);
            }

            //是否接受迭代点
            if (rou > delta_yeta)
            {
                for (int i = 0; i < scope + 2; i++)
                {
                    x[i] = x_Iter[i];
                }
            }


            return x;
        }

        //柯西点法求解子问题的解
        static double[] CalCauchyPoint(double delta, double[] den, double[,] den2)
        {
            double ghg = 0;
            double[] ghg1 = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                double ghg2 = 0;
                for (int j = 0; j < scope + 2; j++)
                {
                    ghg2 += den[j] * den2[j, i];
                }
                ghg1[i] = ghg2;
            }
            for (int i = 0; i < scope + 2; i++)
            {
                ghg += ghg1[i] * den[i];
            }
            double g_Norm = 0;
            for (int i = 0; i < scope + 2; i++)
            {
                g_Norm += den[i] * den[i];
            }
            g_Norm = Math.Sqrt(g_Norm);
            double tao = 0;
            if (ghg <= 0)
            {
                tao = 1;
            }
            else
            {
                tao = Math.Min(1, (Math.Pow(g_Norm, 3))/ (delta * ghg));  
            }

            double[] cauchyPoint = new double[scope + 2];
            for (int i = 0; i < scope + 2; i++)
            {
                cauchyPoint[i] = -tao * delta * den[i] / g_Norm;
            }
            return cauchyPoint;
        }
        #endregion

        //
        #region 共轭梯度法
        static double g0 = 0;
        static double[] d0 = new double[scope + 2];

        //此函数用于迭代共轭梯度法中的x,选用的是Fletcher-Reeves公式
        static double[] CalConjGrads(int scope, double[] x, int l)
        {
            double[] xm = new double[scope + 2];
            double[] g = new double[scope + 2];
            double[] d = new double[scope + 2];
            double[][] G = GetG(scope);
            for (int i = 0; i < scope + 2; i++)
            {
                for (int j = 0; j < scope + 2; j++)
                {
                    g[i] += G[i][j] * x[j];
                    if (l == 0)
                    {
                        d[i] -= G[i][j] * x[j];
                    }
                }
            }
            double p = 0; double q = 0; double[] q1 = new double[scope + 2];
            if (l != 0)
            {
                double b = 0;
                for (int i = 0; i < scope + 2; i++)
                {
                    p += g[i] * g[i];
                }
                b = p / g0; g0 = p;
                for (int i = 0; i < scope + 2; i++)
                {
                    d[i] = -g[i] + b * d0[i];
                }
                for (int i = 0; i < scope + 2; i++)
                {
                    for (int j = 0; j < scope + 2; j++)
                    {
                        q1[i] += d[j] * G[j][i];
                    }
                    q += d[i] * q1[i];
                }
                p /= q;
                for (int i = 0; i < scope + 2; i++)
                {
                    d0[i] = d[i];
                }
                for (int i = 0; i < scope + 2; i++)
                {
                    xm[i] = x[i] + p * d[i];
                }
            }
            else
            {
                for (int i = 0; i < scope + 2; i++)
                {
                    p += g[i] * g[i];
                    for (int j = 0; j < scope + 2; j++)
                    {
                        q1[i] += d[j] * G[j][i];
                    }
                    q += d[i] * q1[i];
                }
                g0 = p;
                p /= q;
                for (int i = 0; i < scope + 2; i++)
                {
                    d0[i] = d[i];
                }
                for (int i = 0; i < scope + 2; i++)
                {
                    xm[i] = x[i] + p * d[i];
                }
            }
            return xm;
        }

        static double[][] GetG(int scope)
        {
            double[][] G = new double[scope + 2][];
            //构造矩阵G
            for (int i = 0; i < scope + 2; i++)
            {
                G[i] = new double[scope + 2];
                for (int j = 0; j < scope + 2; j++)
                {
                    if (i == 0)
                    {
                        G[i][0] = 6; G[i][1] = -2; G[i][2] = -2;
                    }
                    else if (i == 1)
                    {
                        G[i][0] = -2; G[i][1] = 12; G[i][2] = -4; G[i][2] = -2;
                    }
                    else if (i == scope)
                    {
                        G[i][scope - 2] = -2; G[i][scope - 1] = -4; G[i][scope] = 12; G[i][scope + 1] = -2;
                    }
                    else if (i == scope + 1)
                    {
                        G[i][scope - 1] = -2; G[i][scope] = -2; G[i][scope + 1] = 6;
                    }
                    else
                    {
                        G[i][i] = 18; G[i][i - 1] = -4; G[i][i - 2] = -2; G[i][i + 1] = -4; G[i][i + 2] = -2;
                    }
                }
            }
            return G;
        }
        #endregion
    }
}