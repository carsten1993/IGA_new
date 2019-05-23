using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace IGA.Components._01._Model
{
    public class Create_Shell_From_Surface : GH_Component
    {

        public Create_Shell_From_Surface()
          : base("Create Shell From Surface", "create_shell_from_surfacekname",
              "Creates a shell, given surface points, thickness, mesh dimensions and degrees",
              "IGA", "01. Model")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Corner points", "cp", "Corner points of box", GH_ParamAccess.list);
            pManager.AddNumberParameter("Thickness", "t", "Thickness", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Mesh dim", "dim", "Dimensions of mesh [U, V, W]", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Degrees", "deg", "Degrees [p, q, r]", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Patch", "patch", "Patch", GH_ParamAccess.item);
            pManager.AddTextParameter("Info", "info", "Info", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {

            ///////////////////////////////////////////////// INPUT ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            List<Point3d> cp = new List<Point3d>();
            double t = 0;
            List<int> dim = new List<int>();
            List<int> deg = new List<int>();


            if (!DA.GetDataList(0, cp)) return;
            if (!DA.GetData(1, ref t)) return;
            if (!DA.GetDataList(2, dim)) return;
            if (!DA.GetDataList(3, deg)) return;

            int p = deg[0];
            int q = deg[1];
            int r = deg[2];
            int U = dim[0];
            int V = dim[1];
            int W = dim[2];

            if (t <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Thhickness must be greater than zero");
                return;
            }

            /////////////////////////////////////////////// FUNCTIONS //////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            List<double> knotsCsi = CreateKnotvector(p, U);
            List<double> knotsEta = CreateKnotvector(q, V);
            List<double> knotsZeta = CreateKnotvector(r, W);

            int n = knotsCsi.Count - (p + 1);
            int m = knotsEta.Count - (q + 1);
            int l = knotsZeta.Count - (r + 1);

            cp.AddRange(UpdateControlPoints(cp, t, l));
            List<List<List<Point3d>>> B = CreateControlPointList(cp, n, m, l);
            (List<List<int>> INN, List<List<int>> IEN) = CreateINN_IEN(p, q, r, n, m, l);

            Patch patch = new Patch(B, knotsCsi, knotsEta, knotsZeta, p, q, r, INN, IEN);
            string info = patch.GetMeshInfo();

            //////////////////////////////////////////////// OUTPUT ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            DA.SetData(0, patch);
            DA.SetData(1, info);

        }

        //////////////////////////////////////////////// METHODS ///////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        List<double> CreateKnotvector(int deg, int dim)
        {
            List<double> knots = new List<double>();

            for (int i = 0; i <= deg; i++)
            {
                knots.Add(0);
            }

            double dDim = dim;
            double step = 1 / dDim;
            for (int i = 1; i < dim; i++)
            {
                knots.Add(i * step);
            }

            for (int i = 0; i <= deg; i++)
            {
                knots.Add(1);
            }

            return knots;
        }

        List<Point3d> UpdateControlPoints(List<Point3d> controlPoints, double t, int l)
        {
            List<Point3d> ncp = new List<Point3d>();
            double dl = l;
            double step = t / dl;

            for (int i = 1; i < l; i++)
            {
                foreach (Point3d p in controlPoints)
                {
                    ncp.Add(new Point3d(p.X, p.Y, p.Z + i * step));
                }
            }


            return ncp;
        }

        List<List<List<Point3d>>> CreateControlPointList(List<Point3d> controlPoints, int n, int m, int l)
        {
            List<List<List<Point3d>>> B = new List<List<List<Point3d>>>();

            for (int k = 0; k < l; k++)
            {
                List<List<Point3d>> B_i_j = new List<List<Point3d>>();
                for (int j = 0; j < m; j++)
                {
                    List<Point3d> B_i = new List<Point3d>();
                    for (int i = 0; i < n; i++)
                    {
                        B_i.Add(controlPoints[k * n * m + j * n + i]);
                    }
                    B_i_j.Add(B_i);
                }
                B.Add(B_i_j);
            }

            return B;
        }

        (List<List<int>>, List<List<int>>) CreateINN_IEN(int p, int q, int r, int n, int m, int l)
        {
            List<List<int>> INN = new List<List<int>>();
            List<List<int>> IEN = new List<List<int>>();

            int A = 0, B;
            for (int k = 1; k <= l; k++)
            {
                for (int j = 1; j <= m; j++)
                {
                    for (int i = 1; i <= n; i++)
                    {
                        List<int> tempINN = new List<int>();
                        tempINN.Add(i);
                        tempINN.Add(j);
                        tempINN.Add(k);
                        INN.Add(tempINN);
                        A++;

                        if (i >= p + 1 && j >= q + 1 && k >= r + 1)
                        {
                            List<int> tempIEN = new List<int>();
                            for (int kloc = 0; kloc <= r; kloc++)
                            {
                                for (int jloc = 0; jloc <= q; jloc++)
                                {
                                    for (int iloc = 0; iloc <= p; iloc++)
                                    {
                                        B = A - kloc * n * m - jloc * n - iloc;
                                        tempIEN.Add(B);
                                    }
                                }
                            }
                            IEN.Add(tempIEN);
                        }
                    }
                }
            }

            return (INN, IEN);
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
                return IGA.Properties.Resources.surfacemesh;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("fd5fd2eb-44ac-4cb4-b202-d7d0a87468b2"); }
        }
    }
}