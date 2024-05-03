namespace MGroup.Solvers.DDM.Output
{
	using System.Collections.Generic;
	using System.IO;
	using System.Linq;

	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Meshes.Output.VTK;
	using MGroup.Solvers.DDM.Partitioning;

	public class PartitioningWriter
	{
		private readonly string outputDirectory;
		private readonly int dimension;

		public PartitioningWriter(string outputDirectory, int dimension)
		{
			this.outputDirectory = outputDirectory;
			this.dimension = dimension;
		}

		public void PlotPartitioning(IModel model, IPartitioner partitioner, int iteration = int.MinValue)
		{
			List<INode> nodes = model.EnumerateNodes().ToList();
			var allElements = new HashSet<IElementType>();
			foreach (ISubdomain subdomain in model.EnumerateSubdomains())
			{
				allElements.UnionWith(subdomain.EnumerateElements());
			}
			var outputMesh = new VtkMeshDiscontinuous(nodes, allElements.ToList());
			var subdomainIDs = new double[outputMesh.VtkPoints.Count];
			var clusterIDs = new double[outputMesh.VtkPoints.Count];
			for (int e = 0; e < outputMesh.OriginalElements.Count; ++e)
			{
				int subdomainID = partitioner.GetSubdomainOfElement(outputMesh.OriginalElements[e].ID);
				int clusterID = partitioner.GetClusterOfSubdomain(subdomainID);
				VtkCell vtkCell = outputMesh.VtkCells[e];
				foreach (VtkPoint vtkPoint in vtkCell.Vertices)
				{
					subdomainIDs[vtkPoint.ID] = subdomainID;
					clusterIDs[vtkPoint.ID] = clusterID;
				}
			}

			string extension = iteration == int.MinValue ? ".vtk" : $"_{iteration}.vtk";
			string path = Path.Combine(outputDirectory, "partitioning" + extension);

			using (var writer = new VtkFileWriter(path, dimension))
			{
				writer.WriteMesh(outputMesh);
				writer.WriteScalarField("subdomainIDs", subdomainIDs);
				writer.WriteScalarField("clusterIDs", clusterIDs);
			}
		}
	}
}
