namespace MGroup.Solvers.DDM.Tests.ScalabilityAnalysis
{
	using MGroup.Constitutive.Structural;
	using MGroup.Constitutive.Structural.Planar;
	using MGroup.Environments;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.Tests.Commons;

	public class CantilevelBeam2D : IModelBuilder
	{
		public double[] DomainLengthPerAxis { get; set; } = { 8, 2 };

		public int[] NumElementsPerAxis { get; set; }

		public int[] NumSubdomainsPerAxis { get; set; }

		public int[] NumClustersPerAxis { get; set; } = { 2, 2 };

		public int[] SubdomainSizePerElementSize
		{
			get
			{
				var result = new int[2];
				for (int i = 0; i < result.Length; i++)
				{
					if (NumElementsPerAxis[i] % NumSubdomainsPerAxis[i] != 0)
					{
						throw new ArgumentException("Elements per axis must be a multple of subdomains per axis");
					}
					result[i] = NumElementsPerAxis[i] / NumSubdomainsPerAxis[i];
				}
				return result;
			}
		}

		public double YoungModulus { get; set; } = 2.1E7;

		public double EndPointLoad { get; set; } = -2E4;

		public (IModel model, ComputeNodeTopology nodeTopology) CreateMultiSubdomainModel()
		{
			UniformDdmModelBuilder2D builder = CreateDefaultModelBuilder();
			builder.NumSubdomains = NumSubdomainsPerAxis;
			builder.NumClusters = NumClustersPerAxis;
			return builder.BuildMultiSubdomainModel();
		}

		public IModel CreateSingleSubdomainModel()
		{
			UniformDdmModelBuilder2D builder = CreateDefaultModelBuilder();
			builder.NumSubdomains = new int[] { 1, 1 };
			builder.NumClusters = new int[] { 1, 1 };
			return builder.BuildSingleSubdomainModel();
		}

		public ICornerDofSelection GetCornerDofs(IModel model) => UniformDdmModelBuilder2D.FindCornerDofs(model, 2);

		public (List<int[]> numElements, int[] numSubdomains) GetParametricConfigConstNumSubdomains()
		{
			int[] numSubdomains = { 16, 4 };
			var numElements = new List<int[]>();

			numElements.Add(new int[] { 32, 8 });
			numElements.Add(new int[] { 64, 16 });
			numElements.Add(new int[] { 128, 32 });
			numElements.Add(new int[] { 256, 64 });
			//numElements.Add(new int[] { 512, 128 });
			//numElements.Add(new int[] { 1024, 256 });
			//numElements.Add(new int[] { 2048, 512 });

			return (numElements, numSubdomains);
		}

		public (int[] numElements, List<int[]> numSubdomains) GetParametricConfigConstNumElements()
		{
			int[] numElements = { 512, 128 };
			var numSubdomains = new List<int[]>();

			//numSubdomains.Add(new int[] { 4, 1 });
			numSubdomains.Add(new int[] { 8, 2 });
			numSubdomains.Add(new int[] { 16, 4 });
			numSubdomains.Add(new int[] { 32, 8 });
			numSubdomains.Add(new int[] { 64, 16 });
			//numSubdomains.Add(new int[] { 128, 32 });
			//numSubdomains.Add(new int[] { 256, 64 });

			return (numElements, numSubdomains);
		}

		public (List<int[]> numElements, List<int[]> numSubdomains) GetParametricConfigConstSubdomainPerElementSize()
		{
			var numElements = new List<int[]>();
			var numSubdomains = new List<int[]>();

			numElements.Add(new int[] { 32, 8 });
			numElements.Add(new int[] { 64, 16 });
			numElements.Add(new int[] { 128, 32 });
			numElements.Add(new int[] { 256, 64 });
			numElements.Add(new int[] { 512, 128 });
			//numElements.Add(new int[] { 1024, 256 });
			//numElements.Add(new int[] { 2048, 512 });

			numSubdomains.Add(new int[] { 4, 1 });
			numSubdomains.Add(new int[] { 8, 2 });
			numSubdomains.Add(new int[] { 16, 4 });
			numSubdomains.Add(new int[] { 32, 8 });
			numSubdomains.Add(new int[] { 64, 16 });
			numSubdomains.Add(new int[] { 128, 32 });
			//numSubdomains.Add(new int[] { 256, 64 });

			return (numElements, numSubdomains);
		}

		private UniformDdmModelBuilder2D CreateDefaultModelBuilder()
		{
			var builder = new UniformDdmModelBuilder2D();
			builder.MinCoords = new double[2];
			builder.MaxCoords = DomainLengthPerAxis;
			builder.NumElementsTotal = NumElementsPerAxis;

			builder.MaterialHomogeneous = new ElasticMaterial2D(YoungModulus, 0.3, StressState2D.PlaneStress);
			builder.Thickness = 1.0;

			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
			builder.DistributeLoadAtNodes(
				UniformDdmModelBuilder2D.BoundaryRegion.RightSide, StructuralDof.TranslationY, EndPointLoad);

			return builder;
		}
	}
}
