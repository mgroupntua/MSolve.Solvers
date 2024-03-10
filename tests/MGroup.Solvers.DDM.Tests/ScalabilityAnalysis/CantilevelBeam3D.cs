using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Environments;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.Tests.Commons;

namespace MGroup.Solvers.DDM.Tests.ScalabilityAnalysis
{
	public class CantilevelBeam3D : IModelBuilder
	{
		public double[] DomainLengthPerAxis { get; set; } = { 6, 2, 2 };

		public int[] NumElementsPerAxis { get; set; }

		public int[] NumSubdomainsPerAxis { get; set; }

		public int[] NumClustersPerAxis { get; set; } = { 2, 2, 1 };

		public int[] SubdomainSizePerElementSize
		{
			get
			{
				var result = new int[3];
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
			UniformDdmModelBuilder3D builder = CreateDefaultModelBuilder();
			builder.NumSubdomains = NumSubdomainsPerAxis;
			builder.NumClusters = NumClustersPerAxis;
			return builder.BuildMultiSubdomainModel();
		}

		public IModel CreateSingleSubdomainModel()
		{
			UniformDdmModelBuilder3D builder = CreateDefaultModelBuilder();
			builder.NumSubdomains = new int[] { 1, 1 };
			builder.NumClusters = new int[] { 1, 1 };
			return builder.BuildSingleSubdomainModel();
		}

		public ICornerDofSelection GetCornerDofs(IModel model) => UniformDdmModelBuilder3D.FindCornerDofs(model);

		public (List<int[]> numElements, int[] numSubdomains) GetParametricConfigConstNumSubdomains()
		{
			int[] numSubdomains = { 12, 4, 4 };
			var numElements = new List<int[]>();

			numElements.Add(new int[] { 24, 8, 8 });
			numElements.Add(new int[] { 36, 12, 12 });
			numElements.Add(new int[] { 72, 24, 24 });
			//numElements.Add(new int[] { 144, 48, 48 });
			//numElements.Add(new int[] { 288, 96, 96 });
			//numElements.Add(new int[] { 576, 192, 192 });
			//numElements.Add(new int[] { 864, 288, 288 });

			return (numElements, numSubdomains);
		}

		public (int[] numElements, List<int[]> numSubdomains) GetParametricConfigConstNumElements()
		{
			int[] numElements = { 576, 192, 192 };
			var numSubdomains = new List<int[]>();

			numSubdomains.Add(new int[] { 3, 1, 1 });
			numSubdomains.Add(new int[] { 6, 2, 2 });
			numSubdomains.Add(new int[] { 9, 3, 3 });
			numSubdomains.Add(new int[] { 12, 4, 4 });
			numSubdomains.Add(new int[] { 18, 6, 6 });
			numSubdomains.Add(new int[] { 24, 8, 8 });
			numSubdomains.Add(new int[] { 36, 12, 12 });

			return (numElements, numSubdomains);
		}

		public (List<int[]> numElements, List<int[]> numSubdomains) GetParametricConfigConstSubdomainPerElementSize()
		{
			var numElements = new List<int[]>();
			var numSubdomains = new List<int[]>();

			numElements.Add(new int[] { 27, 9, 9 });
			numElements.Add(new int[] { 54, 18, 18 });
			numElements.Add(new int[] { 81, 27, 27 });
			numElements.Add(new int[] { 108, 36, 36 });
			numElements.Add(new int[] { 162, 54, 54 });
			numElements.Add(new int[] { 216, 72, 72 });
			numElements.Add(new int[] { 324, 108, 108 });

			numSubdomains.Add(new int[] { 3, 1, 1 });
			numSubdomains.Add(new int[] { 6, 2, 2 });
			numSubdomains.Add(new int[] { 9, 3, 3 });
			numSubdomains.Add(new int[] { 12, 4, 4 });
			numSubdomains.Add(new int[] { 18, 6, 6 });
			numSubdomains.Add(new int[] { 24, 8, 8 });
			numSubdomains.Add(new int[] { 36, 12, 12 });

			return (numElements, numSubdomains);
		}

		private UniformDdmModelBuilder3D CreateDefaultModelBuilder()
		{
			var builder = new UniformDdmModelBuilder3D();
			builder.MinCoords = new double[3];
			builder.MaxCoords = DomainLengthPerAxis;
			builder.NumElementsTotal = NumElementsPerAxis;

			builder.MaterialHomogeneous = new ElasticMaterial3D(YoungModulus, 0.3);

			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinX, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinX, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinX, StructuralDof.TranslationZ, 0.0);
			builder.DistributeLoadAtNodes(
				UniformDdmModelBuilder3D.BoundaryRegion.MaxX, StructuralDof.TranslationY, EndPointLoad);

			return builder;
		}
	}
}
