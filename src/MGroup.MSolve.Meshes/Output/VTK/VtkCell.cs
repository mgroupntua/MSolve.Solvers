using System;
using System.Collections.Generic;
using MGroup.MSolve.Discretization;

namespace MGroup.MSolve.Meshes.Output.VTK
{
    /// <summary>
    /// Cell used to represent VTK grids.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class VtkCell
    {
        public VtkCell(CellType cellType, IReadOnlyList<VtkPoint> vertices)
        {
            try
            {
                int code = VtkCellTypes.GetVtkCellCode(cellType);
                this.Code = code;
                this.Vertices = vertices;
            } 
            catch (KeyNotFoundException)
            {
                throw new NotImplementedException("Cannot plot elements of type " + cellType);
            }
        }

        public int Code { get; }

        public IReadOnlyList<VtkPoint> Vertices { get; }
    }
}
