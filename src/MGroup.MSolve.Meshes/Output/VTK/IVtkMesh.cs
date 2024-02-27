using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.MSolve.Meshes.Output.VTK
{
    public interface IVtkMesh
    {
        IReadOnlyList<VtkCell> VtkCells { get; } //TODO: Make this a dictionary

        IReadOnlyList<VtkPoint> VtkPoints { get; } //TODO: Make this a dictionary
    }
}
