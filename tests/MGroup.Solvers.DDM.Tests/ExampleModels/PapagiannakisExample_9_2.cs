using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.Planar;
using MGroup.Environments;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.Tests.Commons;
using MGroup.Solvers.Results;

// Subdomains:
// /|                 |\
// /||-------|-------||\
// /||  (2)  |  (3)  ||\
// /||   E1  |   E0  ||\
// /||-------|-------||\
// /||  (0)  |  (1)  ||\
// /||   E1  |   E0  ||\
// /||-------|-------||\
// /|				  |\
namespace MGroup.Solvers.DDM.Tests.ExampleModels
{
	/// <summary>
	/// See Papagiannakis Bachelor thesis, pages 148-161
	/// </summary>
	public static class PapagiannakisExample_9_2
	{
		private const double E0 = 2.1E7, v = 0.3, thickness = 0.3;
		private const double load = 100;

		public static double[] MinCoords => new double[] { 0, 0 };

		public static double[] MaxCoords => new double[] { 3, 1.5 };

		public static int[] NumElements => new int[] { 20, 20 };

		public static int[] NumSubdomains => new int[] { 2, 2 };

		public static int[] NumClusters => new int[] { 1, 1 };

		public static int NumTotalDofs => 882;

		public static Model CreateSingleSubdomainModel(double stiffnessRatio)
		{
			UniformDdmModelBuilder2D builder = CreateDefaultModelBuilder(stiffnessRatio);
			builder.NumSubdomains = new int[] { 1, 1 };
			builder.NumClusters = new int[] { 1, 1 };
			return builder.BuildSingleSubdomainModel();
		}

		public static (Model model, ComputeNodeTopology topology) CreateMultiSubdomainModel(double stiffnessRatio)
		{
			UniformDdmModelBuilder2D builder = CreateDefaultModelBuilder(stiffnessRatio);
			builder.NumSubdomains = NumSubdomains;
			builder.NumClusters = NumClusters;
			return builder.BuildMultiSubdomainModel();
		}

		public static ICornerDofSelection GetCornerDofs(IModel model) => UniformDdmModelBuilder2D.FindCornerDofs(model, 2);

		public static NodalResults SolveWithSkylineSolver(Model model) => PapagiannakisExample_8.SolveWithSkylineSolver(model);

		private static UniformDdmModelBuilder2D CreateDefaultModelBuilder(double stiffnessRatio)
		{
			var builder = new UniformDdmModelBuilder2D();
			builder.MinCoords = MinCoords;
			builder.MaxCoords = MaxCoords;
			builder.NumElementsTotal = NumElements;

			double E1 = E0 / stiffnessRatio;
			builder.GetMaterialPerElementIndex = elementIdx =>
			{
				if (elementIdx[0] < 10)
				{
					return new ElasticMaterial2D(E1, v, StressState2D.PlaneStress);
				}
				else
				{
					return new ElasticMaterial2D(E0, v, StressState2D.PlaneStress);
				}
			};
			builder.Thickness = thickness;

			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.RightSide, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.RightSide, StructuralDof.TranslationY, 0.0);
			builder.DistributeLoadAtNodes(UniformDdmModelBuilder2D.BoundaryRegion.UpperSide, StructuralDof.TranslationY, load);

			return builder;
		}
	}
}
