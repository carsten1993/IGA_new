using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using IGA.Class;

namespace IGA.Components._04._Loads
{
    public class Create_Load : GH_Component
    {

        public Create_Load()
          : base("Create Load", "create_load",
              "Applies a load to surfaces",
              "IGA", "04. Loads")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Patch", "patch", "Patch", GH_ParamAccess.item);
            pManager.AddGenericParameter("Surfaces", "surfaces", "List of surfaces to be loaded", GH_ParamAccess.list);
            pManager.AddGenericParameter("Load vector", "load", "Loading as a vector [X, Y, Z]", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("External load vector", "P1", "External load vector for the given load", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {

            ///////////////////////////////////////////////// INPUT ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            Patch patch = new Patch();
            List<IsoGeoSurface> surfaces = new List<IsoGeoSurface>();
            GH_Vector load = new GH_Vector();

            if (!DA.GetData(0, ref patch)) return;
            if (!DA.GetDataList(1, surfaces)) return;
            if (!DA.GetData(2, ref load)) return;

            if (surfaces.Count <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No nodes are selected");

            }

            (int nel, int nnp, int nen) = patch.GetSize();

            /////////////////////////////////////////////// FUNCTIONS //////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            Vector<double> P1 = DenseVector.OfArray(new double [3 * nnp]);

            foreach (IsoGeoSurface surface in surfaces)
            {
                P1 += CreateExternalLoadVector(patch, surface, load);
            }

            //////////////////////////////////////////////// OUTPUT ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            DA.SetData(0, P1);

        }

        //////////////////////////////////////////////// METHODS ///////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        List<int> Sort(List<int> idLists)
        {
            List<int> id = new List<int>();
            id.Add(idLists[0]);
            int temp1, temp2;
            int c;

            for (int i = 1; i < idLists.Count; i++)
            {
                temp1 = idLists[i];
                c = 0;

                for (int j = 0; j < id.Count; j++)
                {
                    temp2 = id[j];
                    if (temp1 == temp2)
                    {
                        c++;
                    }
                }

                if (c == 0)
                {
                    id.Add(idLists[i]);
                }
            }

            return id;
        }

        (List<double>, List<double>) Gauss(int ngp)
        {

            List<double> list1 = new List<double>(); //Sample Coord in order(exterior-middle-center)
            List<double> list2 = new List<double>(); //Weights


            switch (ngp)
            {
                case 1:
                    list1.Add(0.0);
                    list2.Add(2.0);
                    return (list1, list2);

                case 2:
                    double C = 0.5773502691;
                    list1.Add(-C);
                    list1.Add(C);
                    list2.Add(1.0);
                    list2.Add(1.0);
                    return (list1, list2);

                case 3:
                    double C1 = 0.774596669241483;
                    double C2 = 0.0;
                    double W1 = 0.555555555555555;
                    double W2 = 0.88888888888888;
                    list1.Add(-C1);
                    list1.Add(C2);
                    list1.Add(C1);
                    list2.Add(W1);
                    list2.Add(W2);
                    list2.Add(W1);
                    return (list1, list2);
                case 4:
                    double C11 = 0.861136311594053;
                    double C22 = 0.339981043584856;
                    double W11 = 0.347854845137454;
                    double W22 = 0.652145154862546;
                    list1.Add(-C11);
                    list1.Add(-C22);
                    list1.Add(C22);
                    list1.Add(C11);
                    list2.Add(W11);
                    list2.Add(W22);
                    list2.Add(W22);
                    list2.Add(W11);
                    return (list1, list2);

            }
            return (list1, list2);
        }

        double BasisFunction(List<double> knotsCsi, int i, int p, double csi)
        {
            double sum = 0;

            if (p == 0)
            {
                if (csi >= knotsCsi[i] && csi <= knotsCsi[i + 1])
                {
                    return 1;
                }
                else return 0;
            }
            else
            {
                double d1 = knotsCsi[i + p] - knotsCsi[i];
                double d2 = knotsCsi[i + p + 1] - knotsCsi[i + 1];

                if (d1 != 0 && d2 == 0)
                {
                    double e1 = csi - knotsCsi[i];
                    sum = e1 * BasisFunction(knotsCsi, i, p - 1, csi) / d1;
                }
                else if (d1 == 0 && d2 != 0)
                {
                    double e2 = knotsCsi[i + p + 1] - csi;
                    sum = e2 * BasisFunction(knotsCsi, i + 1, p - 1, csi) / d2;
                }
                else if (d1 != 0 && d2 != 0)
                {
                    double e1 = csi - knotsCsi[i];
                    double e2 = knotsCsi[i + p + 1] - csi;

                    if (csi == knotsCsi[i + 1])
                    {
                        sum = (e1 * BasisFunction(knotsCsi, i, p - 1, csi) / d1);
                    }
                    else
                    {
                        sum = (e1 * BasisFunction(knotsCsi, i, p - 1, csi) / d1) + (e2 * BasisFunction(knotsCsi, i + 1, p - 1, csi) / d2);
                    }
                }
                else return 0;

                return sum;
            }
        }

        double DerivativeBasisFunction(List<double> knotsCsi, int i, int p, double csi)
        {
            double sum = 0;

            if (p == 0)
            {
                if (csi >= knotsCsi[i] && csi <= knotsCsi[i + 1])
                {
                    return 1;
                }
                else return 0;
            }
            else
            {
                double d1 = knotsCsi[i + p] - knotsCsi[i];
                double d2 = knotsCsi[i + p + 1] - knotsCsi[i + 1];

                if (d1 != 0 && d2 == 0)
                {
                    double e1 = p;
                    sum = e1 * BasisFunction(knotsCsi, i, p - 1, csi) / d1;
                }
                else if (d1 == 0 && d2 != 0)
                {
                    double e2 = p;
                    sum = -e2 * BasisFunction(knotsCsi, i + 1, p - 1, csi) / d2;
                }
                else if (d1 != 0 && d2 != 0)
                {
                    double e1 = p;
                    double e2 = p;

                    if (csi == knotsCsi[i + 1])
                    {
                        sum = (e1 * BasisFunction(knotsCsi, i, p - 1, csi) / d1);
                    }
                    else
                    {
                        sum = (e1 * BasisFunction(knotsCsi, i, p - 1, csi) / d1) - (e2 * BasisFunction(knotsCsi, i + 1, p - 1, csi) / d2);
                    }
                }
                else return 0;

                return sum;
            }
        }

        List<int> GetNURBSCoord(List<List<int>> INN, List<List<int>> IEN, int e)
        {
            int ni = INN[IEN[e - 1][0] - 1][0];
            int nj = INN[IEN[e - 1][0] - 1][1];
            int nk = INN[IEN[e - 1][0] - 1][2];

            return new List<int>() { ni, nj, nk };
        }

        List<double> GetParametricCoord(List<int> nCoord, List<double> knotsCsi, List<double> knotsEta, List<double> knotsZeta, List<double> csiTilde)
        {
            int ni = nCoord[0];
            int nj = nCoord[1];
            int nk = nCoord[2];

            double csi = ((knotsCsi[ni] - knotsCsi[ni - 1]) * csiTilde[0] + (knotsCsi[ni] + knotsCsi[ni - 1])) / 2;
            double eta = ((knotsEta[nj] - knotsEta[nj - 1]) * csiTilde[1] + (knotsEta[nj] + knotsEta[nj - 1])) / 2;
            double zeta = ((knotsZeta[nk] - knotsZeta[nk - 1]) * csiTilde[2] + (knotsZeta[nk] + knotsZeta[nk - 1])) / 2;

            return new List<double>() { csi, eta, zeta };
        }
        /*
        (List<List<double>>, List<List<double>>, List<List<double>>) UnivariateBasisFuntion(int p, int q, int r, List<int> nCoord, List<double> pCoord, List<double> knotsCsi, List<double> knotsEta, List<double> knotsZeta)
        {
            List<List<double>> NdN = new List<List<double>>();
            List<double> N = new List<double>();
            List<double> dN = new List<double>();
            List<List<double>> MdM = new List<List<double>>();
            List<double> M = new List<double>();
            List<double> dM = new List<double>();
            List<List<double>> LdL = new List<List<double>>();
            List<double> L = new List<double>();
            List<double> dL = new List<double>();

            int ni = nCoord[0];
            double csi = pCoord[0];
            int nj = nCoord[1];
            double eta = pCoord[1];
            int nk = nCoord[2];
            double zeta = pCoord[2];

            for (int i = ni - p; i <= ni; i++)
            {
                N.Add(BasisFunction(knotsCsi, i - 1, p, csi));
                dN.Add(DerivativeBasisFunction(knotsCsi, i - 1, p, csi));
            }
            NdN.Add(N);
            NdN.Add(dN);

            for (int j = nj - q; j <= nj; j++)
            {
                M.Add(BasisFunction(knotsEta, j - 1, q, eta));
                dM.Add(DerivativeBasisFunction(knotsEta, j - 1, q, eta));
            }
            MdM.Add(M);
            MdM.Add(dM);

            for (int k = nk - r; k <= nk; k++)
            {
                L.Add(BasisFunction(knotsZeta, k - 1, r, zeta));
                dL.Add(DerivativeBasisFunction(knotsZeta, k - 1, r, zeta));
            }
            LdL.Add(L);
            LdL.Add(dL);

            return (NdN, MdM, LdL);
        }
        
        (List<double>, List<List<double>>) TrivariateBasisFunction(int p, int q, int r, List<List<double>> NdN, List<List<double>> MdM, List<List<double>> LdL)
        {
            List<double> R = new List<double>();
            List<List<double>> dR = new List<List<double>>();
            List<double> dR_dCsi = new List<double>();
            List<double> dR_dEta = new List<double>();
            List<double> dR_dZeta = new List<double>();

            for (int k = 0; k <= r; k++)
            {
                for (int j = 0; j <= q; j++)
                {
                    for (int i = 0; i <= p; i++)
                    {
                        R.Add(NdN[0][i] * MdM[0][j] * LdL[0][k]);
                        dR_dCsi.Add(NdN[1][i] * MdM[0][j] * LdL[0][k]);
                        dR_dEta.Add(NdN[0][i] * MdM[1][j] * LdL[0][k]);
                        dR_dZeta.Add(NdN[0][i] * MdM[0][j] * LdL[1][k]);
                    }
                }
            }

            dR.Add(dR_dCsi);
            dR.Add(dR_dEta);
            dR.Add(dR_dZeta);

            return (R, dR);
        }

        (List<double> s_R, List<List<double>> s_dR_dCsi) GetSurfaceTrivariateTrivariateBasisFunction(List<double> R, List<List<double>> dR_dCsi, int e, List<List<int>> IEN, List<int> sf, int type)
        {
            List<double> s_R = new List<double>();
            List<List<double>> s_dR = new List<List<double>>();
            List<double> s_dR_d1 = new List<double>();
            List<double> s_dR_d2 = new List<double>();

            int p1, p2; 
            if (type == 1 || type == 6)
            {
                p1 = 0;
                p2 = 1;
            }
            else if (type == 2 || type == 4)
            {
                p1 = 0;
                p2 = 2;
            }
            else
            {
                p1 = 1;
                p2 = 2;
            }

            for (int i = IEN[e - 1].Count - 1, j = 0; i >= 0; i--, j++)
            {
                foreach (int d in sf)
                {
                    if (IEN[e - 1][i] == d)
                    {
                        s_R.Add(R[j]);
                        s_dR_d1.Add(dR_dCsi[p1][j]);
                        s_dR_d2.Add(dR_dCsi[p2][j]);
                    }
                }
            }

            s_dR.Add(s_dR_d1);
            s_dR.Add(s_dR_d2);

            return (s_R, s_dR);
        }
        */
        List<List<Point3d>> GetSurfaceControlPoints(List<List<List<Point3d>>> B, int type)
        {
            List<List<Point3d>> s_B = new List<List<Point3d>>();
            int xDim = B[0][0].Count;
            int yDim = B[0].Count;
            int zDim = B.Count;

            if (type == 1)
            {
                for (int j = 0; j < yDim; j++)
                {
                    List<Point3d> list1 = new List<Point3d>();
                    for (int i = 0; i < xDim; i++)
                    {
                        Point3d tempPoint = new Point3d(B[0][j][i]);
                        list1.Add(tempPoint);
                    }
                    s_B.Add(list1);
                }
            }
            else if (type == 2)
            {
                for (int k = 0; k < zDim; k++)
                {
                    List<Point3d> list1 = new List<Point3d>();
                    for (int i = 0; i < xDim; i++)
                    {
                        Point3d tempPoint = new Point3d(B[k][0][i]);
                        list1.Add(tempPoint);
                    }
                    s_B.Add(list1);
                }
            }
            else if (type == 3)
            {
                for (int k = 0; k < zDim; k++)
                {
                    List<Point3d> list1 = new List<Point3d>();
                    for (int j = 0; j < yDim; j++)
                    {
                        Point3d tempPoint = new Point3d(B[k][j][xDim - 1]);
                        list1.Add(tempPoint);
                    }
                    s_B.Add(list1);
                }
            }
            else if (type == 4)
            {
                for (int k = 0; k < zDim; k++)
                {
                    List<Point3d> list1 = new List<Point3d>();
                    for (int i = 0; i < xDim; i++)
                    {
                        Point3d tempPoint = new Point3d(B[k][yDim - 1][i]);
                        list1.Add(tempPoint);
                    }
                    s_B.Add(list1);
                }
            }
            else if (type == 5)
            {
                for (int k = 0; k < zDim; k++)
                {
                    List<Point3d> list1 = new List<Point3d>();
                    for (int j = 0; j < yDim; j++)
                    {
                        Point3d tempPoint = new Point3d(B[k][j][0]);
                        list1.Add(tempPoint);
                    }
                    s_B.Add(list1);
                }
            }
            else
            {
                for (int j = 0; j < yDim; j++)
                {
                    List<Point3d> list1 = new List<Point3d>();
                    for (int i = 0; i < xDim; i++)
                    {
                        Point3d tempPoint = new Point3d(B[zDim - 1][j][i]);
                        list1.Add(tempPoint);
                    }
                    s_B.Add(list1);
                }
            }

            return s_B;
        }

        (Matrix<double>, Matrix<double>) CreateGradients(int p1, int p2, List<int> nCoord, List<List<Point3d>> s_B, List<List<double>> dR_dCsi, List<double> knots1, List<double> knots2, int type)
        {
            Matrix<double> dx_dCsi = DenseMatrix.OfArray(new double[2, 2]);
            Matrix<double> dCsi_dCsiTilde = DenseMatrix.OfArray(new double[2, 2]);

            int n1 = nCoord[0];
            int n2 = nCoord[1];

            dCsi_dCsiTilde[0, 0] = (knots1[n1] - knots1[n1 - 1]) / 2;
            dCsi_dCsiTilde[1, 1] = (knots2[n2] - knots2[n2 - 1]) / 2;

            Point3d point;
            int loc = 0;

            if (type == 1 || type == 6)
            {
                for (int j = 0; j <= p2; j++)
                {
                    for (int i = 0; i <= p1; i++)
                    {
                        point = s_B[n2 - p2 + j - 1][n1 - p1 + i - 1];
                        List<double> pointList = new List<double>() { point.X, point.Y };

                        for (int a = 0; a < 2; a++)
                        {
                            for (int b = 0; b < 2; b++)
                            {
                                dx_dCsi[a, b] += dR_dCsi[b][loc] * pointList[a];
                            }
                        }
                        loc++;
                    }
                }
            }
            else if (type == 2 || type == 4)
            {
                for (int j = 0; j <= p2; j++)
                {
                    for (int i = 0; i <= p1; i++)
                    {
                        point = s_B[n2 - p2 + j - 1][n1 - p1 + i - 1];
                        List<double> pointList = new List<double>() { point.X, point.Z };

                        for (int a = 0; a < 2; a++)
                        {
                            for (int b = 0; b < 2; b++)
                            {
                                dx_dCsi[a, b] += dR_dCsi[b][loc] * pointList[a];
                            }
                        }
                        loc++;
                    }
                }
            }
            else
            {
                for (int j = 0; j <= p2; j++)
                {
                    for (int i = 0; i <= p1; i++)
                    {
                        point = s_B[n2 - p2 + j - 1][n1 - p1 + i - 1];
                        List<double> pointList = new List<double>() { point.Y, point.Z };

                        for (int a = 0; a < 2; a++)
                        {
                            for (int b = 0; b < 2; b++)
                            {
                                dx_dCsi[a, b] += dR_dCsi[b][loc] * pointList[a];
                            }
                        }
                        loc++;
                    }
                }
            }



                return (dx_dCsi, dCsi_dCsiTilde);
        }

        (List<double>, List<List<double>>) SurfaceTrivariateBasisFuntion(int p1, int p2, int p3, List<int> n1Coord, List<double> pCoord, List<double> knots1, List<double> knots2, List<double> knots3)
        {
            List<double> R = new List<double>();
            List<List<double>> dR = new List<List<double>>();
            List<double> dR1 = new List<double>();
            List<double> dR2 = new List<double>();

            int n1 = n1Coord[0];
            double csi1 = pCoord[0];
            int n2 = n1Coord[1];
            double csi2 = pCoord[1];
            int n3 = n1Coord[2];
            double csi3 = pCoord[2];
            double bf3 = BasisFunction(knots3, n3 - 1, p3, csi3);

            for (int j = n2 - p2; j <= n2; j++)
            {
                for (int i = n1 - p1; i <= n1; i++)
                {
                    R.Add(BasisFunction(knots1, i - 1, p1, csi1) * BasisFunction(knots2, j - 1, p2, csi2)*bf3);
                    dR1.Add(DerivativeBasisFunction(knots1, i - 1, p1, csi1) * BasisFunction(knots2, j - 1, p2, csi2) * bf3);
                    dR2.Add(BasisFunction(knots1, i - 1, p1, csi1) * DerivativeBasisFunction(knots2, j - 1, p2, csi2) * bf3);
                }
            }
            dR.Add(dR1);
            dR.Add(dR2);

            return (R, dR);
        }

        (List<double>, double) GetShapeFunctions(int p1, int p2, int p3, List<int> n1Coord, List<double> p1Coord, List<double> knots1, List<double> knots2, List<double> knots3, List<List<Point3d>> s_B, int type)
        {
            (List<double> s_R, List<List<double>> s_dR_dCsi) = SurfaceTrivariateBasisFuntion(p1, p2, p3, n1Coord, p1Coord, knots1, knots2, knots3);
            (Matrix<double> s_dx_dCsi, Matrix<double> s_dCsi_dCsiTilde) = CreateGradients(p1, p2, n1Coord, s_B, s_dR_dCsi, knots1, knots2, type);
            Matrix<double> s_J_mat = s_dx_dCsi * s_dCsi_dCsiTilde;
            double s_J = s_J_mat.Determinant();

            return (s_R, s_J);
        }

        Vector<double> CreateGaussPointElementLoadVector(List<double> R, int s_nen, GH_Vector load, double modJ)
        {
            Vector<double> Pgpe = DenseVector.OfArray(new double[3 * s_nen]);

            for (int i = 0; i < s_nen; i++)
            {
                Pgpe[i * 3] = load.Value.X * R[i] * modJ;
                Pgpe[i * 3 + 1] = load.Value.Y * R[i] * modJ;
                Pgpe[i * 3 + 2] = load.Value.Z * R[i] * modJ;
            }

            return Pgpe;
        }

        Vector<double> CreateElementExternalLoadVector(int p, int q, int r, int n, int m, int l, int s_nen, List<int> nCoord, List<double> knotsCsi, List<double> knotsEta, List<double> knotsZeta, List<List<Point3d>> s_B, GH_Vector load, int type)
        {
            Vector<double> P1e = DenseVector.OfArray(new double[3 * s_nen]);

            int p1, p2, p3;
            List<double> knots1, knots2, knots3;
            List<int> n1Coord;
            if (type == 1)
            {
                p1 = p;                
                p2 = q;
                p3 = r;
                knots1 = knotsCsi;
                knots2 = knotsEta;
                knots3 = knotsZeta;
                n1Coord = new List<int>() { nCoord[0], nCoord[1], 1 };
            }
            else if (type == 2)
            {
                p1 = p;
                p2 = r;
                p3 = q;
                knots1 = knotsCsi;
                knots2 = knotsZeta;
                knots3 = knotsEta;
                n1Coord = new List<int>() { nCoord[0], 1, nCoord[2] };
            }
            else if (type == 3)
            {
                p1 = q;
                p2 = r;
                p3 = p;
                knots1 = knotsEta;
                knots2 = knotsZeta;
                knots3 = knotsCsi;
                n1Coord = new List<int>() { n, nCoord[1], nCoord[2] };
            }
            else if (type == 4)
            {
                p1 = p;
                p2 = r;
                p3 = q;
                knots1 = knotsCsi;
                knots2 = knotsZeta;
                knots3 = knotsEta;
                n1Coord = new List<int>() { nCoord[0], m, nCoord[2] };
            }
            else if (type == 5)
            {
                p1 = q;
                p2 = r;
                p3 = p;
                knots1 = knotsEta;
                knots2 = knotsZeta;
                knots3 = knotsCsi;
                n1Coord = new List<int>() { 1, nCoord[1], nCoord[2] };
            }
            else
            {
                p1 = p;
                p2 = q;
                p3 = r;
                knots1 = knotsCsi;
                knots2 = knotsEta;
                knots3 = knotsZeta;
                n1Coord = new List<int>() { nCoord[0], nCoord[1], l };
            }

            int ngp1 = p1 + 1;
            int ngp2 = p2 + 1;
            (List<double> p1TildeCoord, List<double> weight1) = Gauss(ngp1);
            (List<double> p2TildeCoord, List<double> weight2) = Gauss(ngp2);

            for (int j = 0; j < ngp2; j++)
            {
                for (int i = 0; i < ngp1; i++)
                {
                    List<double> csiTilde = new List<double>();
                    if (type == 1 )
                    {
                        csiTilde.Add(p1TildeCoord[i]);
                        csiTilde.Add(p2TildeCoord[j]);
                        csiTilde.Add(-1 );
                    }
                    else if (type == 2 )
                    {
                        csiTilde.Add(p1TildeCoord[i]);
                        csiTilde.Add(-1);
                        csiTilde.Add(p2TildeCoord[j]);
                    }
                    else if (type == 3)
                    {
                        csiTilde.Add(1);
                        csiTilde.Add(p1TildeCoord[i]);
                        csiTilde.Add(p2TildeCoord[j]);
                    }
                    else if (type == 4)
                    {
                        csiTilde.Add(p1TildeCoord[i]);
                        csiTilde.Add(1);
                        csiTilde.Add(p2TildeCoord[j]);
                    }
                    else if (type == 5)
                    {
                        csiTilde.Add(-1);
                        csiTilde.Add(p1TildeCoord[i]);
                        csiTilde.Add(p2TildeCoord[j]);
                    }
                    else
                    {
                        csiTilde.Add(p1TildeCoord[i]);
                        csiTilde.Add(p2TildeCoord[j]);
                        csiTilde.Add(1);
                    }

                    
                    List<double> p1Coord = GetParametricCoord(nCoord, knotsCsi, knotsEta, knotsZeta, csiTilde);
                    (List<double> R, double J) = GetShapeFunctions(p1, p2, p3, n1Coord, p1Coord, knots1, knots2, knots3, s_B, type);
                    double modJ = J * weight1[i] * weight2[j];
                    P1e += CreateGaussPointElementLoadVector(R, s_nen, load, modJ);
                }
            }           

            return P1e;
        }

        Vector<double> CreateExternalLoadVector(Patch patch, IsoGeoSurface surface, GH_Vector load)
        {
            List<List<List<Point3d>>> B = patch.GetControlPoints();
            (List<double> knotsCsi, List<double> knotsEta, List<double> knotsZeta) = patch.GetKnotVectors();
            (int p, int q, int r) = patch.GetDegrees();
            (int n, int m, int l) = patch.GetNML();
            (int nel, int nnp, int nen) = patch.GetSize();
            (List<List<int>> INN, List<List<int>> IEN) = patch.GetINNIEN();

            List<int> sf = surface.GetShapeFunctions();
            int s_nnp = sf.Count;
            List<int> elemList = surface.GetElements();
            List<List<int>> s_IEN = surface.GetIEN();
            int s_nen = s_IEN[0].Count;
            int type = surface.GetSurface();

            List<List<Point3d>> s_B = GetSurfaceControlPoints(B, type);

            Vector<double> P1 = DenseVector.OfArray(new double[3 * nnp]);
            Vector<double> P1e = DenseVector.OfArray(new double[3 * s_nen]);

            int elem = 0;
            foreach (int e in elemList)
            {
                List<int> nCoord = GetNURBSCoord(INN, IEN, e);
                
                int ni = nCoord[0];
                int nj = nCoord[1];
                int nk = nCoord[2];

                
                if (knotsCsi[ni] == knotsCsi[ni - 1] || knotsEta[nj] == knotsEta[nj - 1] || knotsZeta[nk] == knotsZeta[nk - 1])
                {
                    continue;
                }

                P1e = CreateElementExternalLoadVector(p, q, r, n, m, l, s_nen, nCoord, knotsCsi, knotsEta, knotsZeta, s_B, load, type);

                int localI1, localI2, localI3, globalI1, globalI2, globalI3, globalFI;
                for (int i = 0; i < s_nen; i ++)
                {
                    globalFI = s_nen - 1 - i; 

                    localI1 = i*3;
                    localI2 = i*3 + 1;
                    localI3 = i*3 + 2;

                    globalI1 = (s_IEN[elem][globalFI] - 1) * 3;
                    globalI2 = (s_IEN[elem][globalFI] - 1) * 3 + 1;
                    globalI3 = (s_IEN[elem][globalFI] - 1) * 3 + 2;

                    P1[globalI1] += P1e[localI1];
                    P1[globalI2] += P1e[localI2];
                    P1[globalI3] += P1e[localI3];                    
                }

                elem++;
            }

            return P1;
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return IGA.Properties.Resources.pointload;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("0cff849a-eaaa-4cf8-9c48-7314663a4f1a"); }
        }
    }
}