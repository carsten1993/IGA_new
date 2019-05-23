using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

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
            pManager.AddIntegerParameter("Shape functions/control points", "idLists", "List of shape functions/control points to be restrained", GH_ParamAccess.list);
            pManager.AddGenericParameter("Load vector", "load", "Loading as a vector [X, Y, Z]", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {

            ///////////////////////////////////////////////// INPUT ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            List<int> idLists = new List<int>();
            GH_Vector load = new GH_Vector();

            if (!DA.GetDataList(0, idLists)) return;
            if (!DA.GetData(1, ref load)) return;

            if (idLists.Count <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No nodes are selected");

            }

            /////////////////////////////////////////////// FUNCTIONS //////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            List<int> id = Sort(idLists);
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


        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
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