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
        private int type;

        public IsoGeoSurface() { }

        public IsoGeoSurface(List<int> _shapeFunctions, List<int> _elements, int _type)
        {
            shapeFunctions = _shapeFunctions;
            elements = _elements;
            type = _type;
        }

        public List<int> GetShapeFunctions() { return shapeFunctions; }
        public List<int> GetElements() { return elements; }
        public int GetSurface() { return type; }
    }
}
