using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Geometry.Coordinates;

namespace MGroup.MSolve.Meshes.Output.VTK
{
    public class VtkPointWriter : IDisposable
    {
        public const string vtkReaderVersion = "4.1";
        private readonly StreamWriter writer;

        public VtkPointWriter(string filePath)
        {
            this.writer = new StreamWriter(filePath);
            writer.Write("# vtk DataFile Version ");
            writer.WriteLine(vtkReaderVersion);
            writer.WriteLine(filePath);
            writer.Write("ASCII\n\n");
        }

        public void Dispose()
        {
            if (writer != null) writer.Dispose();
        }

        public void WritePoints(IEnumerable<CartesianPoint> points)
        {
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {points.Count()} double");
            foreach (CartesianPoint point in points)
            {
                writer.WriteLine($"{point.X} {point.Y} {point.Z}");
            }
        }

        public void WritePoints(IEnumerable<VtkPoint> points)
        {
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {points.Count()} double");
            foreach (VtkPoint point in points)
            {
                writer.WriteLine($"{point.X} {point.Y} {point.Z}");
            }
        }

        public void WriteScalarField(string fieldName, IReadOnlyDictionary<CartesianPoint, double> pointValues)
        {
            // Points
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {pointValues.Count} double");
            foreach (CartesianPoint point in pointValues.Keys)
            {
                writer.WriteLine($"{point.X} {point.Y} {point.Z}");
            }

            // Values
            writer.Write("\n\n");
            writer.WriteLine($"POINT_DATA {pointValues.Count}");
            writer.WriteLine($"SCALARS {fieldName} double 1");
            writer.WriteLine("LOOKUP_TABLE default");
            foreach (double value in pointValues.Values)
            {
                writer.WriteLine(value);
            }
            writer.WriteLine();
        }

        //TODO: Avoid the duplication.
        public void WriteScalarField(string fieldName, IReadOnlyDictionary<INode, double> nodalValues)
        {
            // Points
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {nodalValues.Count} double");
            foreach (var point in nodalValues.Keys)
            {
                writer.WriteLine($"{point.X} {point.Y} {point.Z}");
            }

            // Values
            writer.Write("\n\n");
            writer.WriteLine($"POINT_DATA {nodalValues.Count}");
            writer.WriteLine($"SCALARS {fieldName} double 1");
            writer.WriteLine("LOOKUP_TABLE default");
            foreach (var value in nodalValues.Values)
            {
                writer.WriteLine(value);
            }
            writer.WriteLine();
        }

        public void WriteVectorField(string fieldName, IReadOnlyDictionary<CartesianPoint, double[]> pointVectors, int dimension)
        {
            // Points
            writer.WriteLine("DATASET UNSTRUCTURED_GRID");
            writer.WriteLine($"POINTS {pointVectors.Count} double");
            foreach (CartesianPoint point in pointVectors.Keys)
            {
                writer.WriteLine($"{point.X} {point.Y} {point.Z}");
            }

            // Values
            writer.Write("\n\n");
            writer.WriteLine($"POINT_DATA {pointVectors.Count}");
            writer.WriteLine($"VECTORS {fieldName} double");
            if (dimension == 2)
            {
                foreach (double[] vector in pointVectors.Values)
                {
                    writer.WriteLine($"{vector[0]} {vector[1]} 0.0");
                }
            }
            else
            {
                Debug.Assert(dimension == 3);
                foreach (double[] vector in pointVectors.Values)
                {
                    writer.WriteLine($"{vector[0]} {vector[1]} {vector[2]}");
                }
            }
            writer.WriteLine();
        }
    }
}
