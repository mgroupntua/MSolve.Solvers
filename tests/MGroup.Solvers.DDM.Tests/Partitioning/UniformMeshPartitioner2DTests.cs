namespace MGroup.Solvers.DDM.Tests.Partitioning
{
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Discretization.Meshes.Structured;
	using MGroup.Solvers.DDM.Output;
	using MGroup.Solvers.DDM.Partitioning;
	using MGroup.Solvers.DDM.Tests.ExampleModels;

	using Xunit;

	[Collection("Sequential")]
	public static class UniformMeshPartitioner2DTests
	{
		//[Fact]
		public static void PlotPartitioning()
		{
			IModel model = Plane2DExample.CreateSingleSubdomainModel();
			model.ConnectDataStructures();

			var mesh = new UniformCartesianMesh2D.Builder(
				Plane2DExample.MinCoords, Plane2DExample.MaxCoords, Plane2DExample.NumElements).BuildMesh();
			var partitioner = new UniformMeshPartitioner2D(mesh, Plane2DExample.NumSubdomains, Plane2DExample.NumClusters);
			partitioner.Partition(model);

			string outputDirectory = @"C:\Users\Serafeim\Desktop\PFETIDP\partitioning\plane2D";
			var writer = new PartitioningWriter(outputDirectory, 2);
			writer.PlotPartitioning(model, partitioner);
		}

		[Fact]
		public static void TestMeshPartitioning()
		{
			var mesh = new UniformCartesianMesh2D.Builder(Plane2DExample.MinCoords, Plane2DExample.MaxCoords, 
				Plane2DExample.NumElements).SetMajorAxis(0).BuildMesh();
			var partitioner = new UniformMeshPartitioner2D(mesh, Plane2DExample.NumSubdomains, Plane2DExample.NumClusters);
			partitioner.Partition(null);

			Dictionary<int, int> expectedPartitioning = Plane2DExample.GetSubdomainsOfElements();
			foreach (int elementID in expectedPartitioning.Keys)
			{
				Assert.Equal(expectedPartitioning[elementID], partitioner.GetSubdomainOfElement(elementID));
			}
		}

		[Fact]
		public static void TestSubdomainClustering()
		{
			var mesh = new UniformCartesianMesh2D.Builder(Plane2DExample.MinCoords, Plane2DExample.MaxCoords, 
				Plane2DExample.NumElements).SetMajorAxis(0).BuildMesh();
			var partitioner = new UniformMeshPartitioner2D(mesh, Plane2DExample.NumSubdomains, Plane2DExample.NumClusters);
			partitioner.Partition(null);

			Dictionary<int, int> expectedClustering = Plane2DExample.GetSubdomainClusters();
			foreach (int subdomainID in expectedClustering.Keys)
			{
				Assert.Equal(expectedClustering[subdomainID], partitioner.GetClusterOfSubdomain(subdomainID));
			}
		}

		[Fact]
		public static void TestSubdomainNeighbors()
		{
			var mesh = new UniformCartesianMesh2D.Builder(Plane2DExample.MinCoords, Plane2DExample.MaxCoords, 
				Plane2DExample.NumElements).SetMajorAxis(0).BuildMesh();
			var partitioner = new UniformMeshPartitioner2D(mesh, Plane2DExample.NumSubdomains, Plane2DExample.NumClusters);
			partitioner.Partition(null);

			Dictionary<int, int[]> expectedSubdomainNeighbors = Plane2DExample.GetSubdomainNeighbors();
			foreach (int subdomainID in expectedSubdomainNeighbors.Keys)
			{
				int[] expectedNeighbors = expectedSubdomainNeighbors[subdomainID].OrderBy(n => n).ToArray();
				int[] computedNeighbors = partitioner.GetNeighboringSubdomains(subdomainID).OrderBy(n => n).ToArray();
				Assert.True(Utilities.AreEqual(expectedNeighbors, computedNeighbors));
			}
		}
	}
}
