using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IGA.Class
{
    class IsoGeoSurface
    {
        private List<int> shapeFunctions, elements;
        private List<List<int>> s_IEN;
        private int type;

        public IsoGeoSurface() { }

        public IsoGeoSurface(List<int> _shapeFunctions, List<int> _elements, List<List<int>> _s_IEN, int _type)
        {
            shapeFunctions = _shapeFunctions;
            elements = _elements;
            s_IEN = _s_IEN;
            type = _type;
        }

        public List<int> GetShapeFunctions() { return shapeFunctions; }
        public List<int> GetElements() { return elements; }
        public List<List<int>> GetIEN() { return s_IEN; }
        public int GetSurface() { return type; }
    }
}
