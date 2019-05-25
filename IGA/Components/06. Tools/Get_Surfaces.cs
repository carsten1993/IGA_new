using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using IGA.Class;
using Rhino.Geometry;

namespace IGA
{
    public class Get_Surfaces : GH_Component
    {

        public Get_Surfaces()
          : base("Get Surfaces", "get_surfaces",
              "Returns the ID of the shape functions/control points at the surfaces for an eight cornered solid",
              "IGA", "06. Tools")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Patch", "patch", "Patch", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Surface 1", "surface1", "Surface 1", GH_ParamAccess.item);
            pManager.AddGenericParameter("Surface 2", "surface2", "Surface 2", GH_ParamAccess.item);
            pManager.AddGenericParameter("Surface 3", "surface3", "Surface 3", GH_ParamAccess.item);
            pManager.AddGenericParameter("Surface 4", "surface4", "Surface 4", GH_ParamAccess.item);
            pManager.AddGenericParameter("Surface 5", "surface5", "Surface 5", GH_ParamAccess.item);
            pManager.AddGenericParameter("Surface 6", "surface6", "Surface 6", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {

            ///////////////////////////////////////////////// INPUT ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            Patch patch = new Patch();

            if (!DA.GetData(0, ref patch)) return;

            (int n, int m, int l) = patch.GetNML();
            (int p, int q, int r) = patch.GetDegrees();
            (List<List<int>> INN, List<List<int>> IEN) = patch.GetINNIEN();

            /////////////////////////////////////////////// FUNCTIONS //////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            List<int> surface1 = patch.GetSurface1();
            List<int> elements1 = GetElementsSuface1(n, m, l, p, q, r);
            List<List<int>> s_IEN1 = GetSurfaceIEN(surface1, elements1, IEN);
            IsoGeoSurface isoGeoSurface1 = new IsoGeoSurface(surface1, elements1, s_IEN1, 1);

            List<int> surface2 = patch.GetSurface2();
            List<int> elements2 = GetElementsSuface2(n, m, l, p, q, r);
            List<List<int>> s_IEN2 = GetSurfaceIEN(surface2, elements2, IEN);
            IsoGeoSurface isoGeoSurface2 = new IsoGeoSurface(surface2, elements2, s_IEN2, 2);

            List<int> surface3 = patch.GetSurface3();
            List<int> elements3 = GetElementsSuface3(n, m, l, p, q, r);
            List<List<int>> s_IEN3 = GetSurfaceIEN(surface3, elements3, IEN);
            IsoGeoSurface isoGeoSurface3 = new IsoGeoSurface(surface3, elements3, s_IEN3, 3);

            List<int> surface4 = patch.GetSurface4();
            List<int> elements4 = GetElementsSuface4(n, m, l, p, q, r);
            List<List<int>> s_IEN4 = GetSurfaceIEN(surface4, elements4, IEN);
            IsoGeoSurface isoGeoSurface4 = new IsoGeoSurface(surface4, elements4, s_IEN4, 4);

            List<int> surface5 = patch.GetSurface5();
            List<int> elements5 = GetElementsSuface5(n, m, l, p, q, r);
            List<List<int>> s_IEN5 = GetSurfaceIEN(surface5, elements5, IEN);
            IsoGeoSurface isoGeoSurface5 = new IsoGeoSurface(surface5, elements5, s_IEN5, 5);

            List<int> surface6 = patch.GetSurface6();
            List<int> elements6 = GetElementsSuface6(n, m, l, p, q, r);
            List<List<int>> s_IEN6 = GetSurfaceIEN(surface6, elements6, IEN);
            IsoGeoSurface isoGeoSurface6 = new IsoGeoSurface(surface6, elements6, s_IEN6, 6);

            //////////////////////////////////////////////// OUTPUT ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            DA.SetData(0, isoGeoSurface1);
            DA.SetData(1, isoGeoSurface2);
            DA.SetData(2, isoGeoSurface3);
            DA.SetData(3, isoGeoSurface4);
            DA.SetData(4, isoGeoSurface5);
            DA.SetData(5, isoGeoSurface6);
        }

        //////////////////////////////////////////////// METHODS ///////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        List<int> GetElementsSuface1(int n, int m, int l, int p, int q, int r)
        {
            List<int> elements = new List<int>();
            int numElemX = n - p;
            int numElemY = m - q;

            for (int i = 1; i <= numElemX * numElemY; i++)
            {
                elements.Add(i);
            }

            return elements;
        }

        List<int> GetElementsSuface2(int n, int m, int l, int p, int q, int r)
        {
            List<int> elements = new List<int>();
            int numElemX = n - p;
            int numElemY = m - q;
            int numElemZ = l - r;

            for (int k = 0; k < numElemZ; k++)
            {
                for (int i = 1; i <= numElemX; i++)
                {
                    elements.Add(i + k * numElemX * numElemY);
                }
            }

            return elements;
        }

        List<int> GetElementsSuface3(int n, int m, int l, int p, int q, int r)
        {
            List<int> elements = new List<int>();
            int numElemX = n - p;
            int numElemY = m - q;
            int numElemZ = l - r;

            for (int k = 0; k < numElemZ; k++)
            {
                for (int j = 1; j <= numElemY; j++)
                {
                    elements.Add(j * numElemX + k * numElemX * numElemY);
                }
            }

            return elements;
        }

        List<int> GetElementsSuface4(int n, int m, int l, int p, int q, int r)
        {
            List<int> elements = new List<int>();
            int numElemX = n - p;
            int numElemY = m - q;
            int numElemZ = l - r;

            for (int k = 0; k < numElemZ; k++)
            {
                for (int i = 1; i <= numElemX; i++)
                {
                    elements.Add(i  + (numElemY - 1) * numElemX + k * numElemX * numElemY);
                }
            }

            return elements;
        }

        List<int> GetElementsSuface5(int n, int m, int l, int p, int q, int r)
        {
            List<int> elements = new List<int>();
            int numElemX = n - p;
            int numElemY = m - q;
            int numElemZ = l - r;

            for (int k = 0; k < numElemZ; k++)
            {
                for (int j = 1; j <= numElemY; j++)
                {
                    elements.Add(j * numElemX - numElemX + 1 + k * numElemX * numElemY);
                }
            }

            return elements;
        }

        List<int> GetElementsSuface6(int n, int m, int l, int p, int q, int r)
        {
            List<int> elements = new List<int>();
            int numElemX = n - p;
            int numElemY = m - q;
            int numElemZ = l - r;

            for (int i = numElemX * numElemY * (numElemZ - 1) + 1; i <= numElemX * numElemY * numElemZ; i++)
            {
                elements.Add(i);
            }

            return elements;
        }

        List<List<int>> GetSurfaceIEN(List<int> surface, List<int> elements, List<List<int>> IEN)
        {
            List<List<int>> s_IEN = new List<List<int>>();

            foreach (int e in elements)
            {
                List<int> s_IENi = new List<int>();
                for (int i = 0 - 1; i < IEN[e - 1].Count; i++)
                {
                    foreach (int d in surface)
                    {
                        if (IEN[e - 1][i] == d)
                        {
                            s_IENi.Add(d);
                        }
                    }
                }
                s_IEN.Add(s_IENi);
            }

            return s_IEN;
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
                return IGA.Properties.Resources.mesh_surface;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("068b069b-9c51-46db-be05-cc470d311b78"); }
        }
    }
}