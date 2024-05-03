namespace MGroup.Solvers.DDM.Tests.ScalabilityAnalysis
{
	using MGroup.Constitutive.Structural;
	using MGroup.Constitutive.Structural.Planar;
	using MGroup.Environments;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.Tests.Commons;

	public class Rve2D : IModelBuilder
	{
		public double[] DomainLengthPerAxis { get; set; } = { 2, 2 };

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

		public double[] CenterDisplacement { get; set; } = { 0.01, -0.02 };

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
			int[] numSubdomains = { 8, 8 };
			var numElements = new List<int[]>();

			numElements.Add(new int[] { 16, 16 });
			numElements.Add(new int[] { 32, 32 });
			numElements.Add(new int[] { 64, 64 });
			numElements.Add(new int[] { 128, 128 });
			//numElements.Add(new int[] { 256, 256 });
			//numElements.Add(new int[] { 512, 512 });
			//numElements.Add(new int[] { 768, 768 });

			return (numElements, numSubdomains);
		}

		public (int[] numElements, List<int[]> numSubdomains) GetParametricConfigConstNumElements()
		{
			int[] numElements = { 256, 256 };
			var numSubdomains = new List<int[]>();

			numSubdomains.Add(new int[] { 4, 4 });
			numSubdomains.Add(new int[] { 8, 8 });
			numSubdomains.Add(new int[] { 16, 16 });
			numSubdomains.Add(new int[] { 32, 32 });
			numSubdomains.Add(new int[] { 64, 64 });
			numSubdomains.Add(new int[] { 128, 128 });
			//numSubdomains.Add(new int[] { 256, 256 });

			return (numElements, numSubdomains);
		}

		public (List<int[]> numElements, List<int[]> numSubdomains) GetParametricConfigConstSubdomainPerElementSize()
		{
			var numElements = new List<int[]>();
			var numSubdomains = new List<int[]>();

			//numSubdomains.Add(new int[] { 4, 4 });
			//numSubdomains.Add(new int[] { 8, 8 });
			//numSubdomains.Add(new int[] { 16, 16 });
			//numSubdomains.Add(new int[] { 32, 32 });
			//numSubdomains.Add(new int[] { 64, 64 });
			numSubdomains.Add(new int[] { 96, 96 });
			//numSubdomains.Add(new int[] { 128, 128 });
			//numSubdomains.Add(new int[] { 256, 256 });

			//numElements.Add(new int[] { 32, 32 });
			//numElements.Add(new int[] { 64, 64 });
			//numElements.Add(new int[] { 128, 128 });
			//numElements.Add(new int[] { 256, 256 });
			//numElements.Add(new int[] { 512, 512 });
			numElements.Add(new int[] { 768, 768 });
			//numElements.Add(new int[] { 1024, 1024 });
			//numElements.Add(new int[] { 2048, 2048 });

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
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.RightSide, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.RightSide, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.UpperSide, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.UpperSide, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.LowerSide, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.LowerSide, StructuralDof.TranslationY, 0.0);

			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.Center, StructuralDof.TranslationX,
				CenterDisplacement[0]);
			builder.PrescribeDisplacement(UniformDdmModelBuilder2D.BoundaryRegion.Center, StructuralDof.TranslationY,
				CenterDisplacement[1]);

			return builder;
		}
	}
}
