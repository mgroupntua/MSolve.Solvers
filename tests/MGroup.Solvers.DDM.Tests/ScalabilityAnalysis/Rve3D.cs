using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Environments;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.Tests.Commons;

namespace MGroup.Solvers.DDM.Tests.ScalabilityAnalysis
{
	public class Rve3D : IModelBuilder
	{
		public double[] DomainLengthPerAxis { get; set; } = { 2, 2, 2 };

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

		public double[] CenterDisplacement { get; set; } = { 0.01, -0.02, 0.02 };

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
			int[] numSubdomains = { 8, 8, 8 };
			var numElements = new List<int[]>();

			numElements.Add(new int[] { 16, 16, 16 });
			numElements.Add(new int[] { 24, 24, 24 });
			numElements.Add(new int[] { 32, 32, 32 });
			numElements.Add(new int[] { 40, 40, 40 });
			//numElements.Add(new int[] { 48, 48, 48 });
			//numElements.Add(new int[] { 56, 56, 56 });
			//numElements.Add(new int[] { 64, 64, 64 });
			//numElements.Add(new int[] { 72, 72, 72 });
			//numElements.Add(new int[] { 80, 80, 80 });

			return (numElements, numSubdomains);
		}

		public (int[] numElements, List<int[]> numSubdomains) GetParametricConfigConstNumElements()
		{
			int[] numElements = { 64, 64, 64 };
			var numSubdomains = new List<int[]>();

			//numSubdomains.Add(new int[] { 2, 2, 2 });
			numSubdomains.Add(new int[] { 4, 4, 4 });
			numSubdomains.Add(new int[] { 8, 8, 8 });
			numSubdomains.Add(new int[] { 16, 16, 16 });
			numSubdomains.Add(new int[] { 32, 32, 32 });
			//numSubdomains.Add(new int[] { 64, 64, 64 });
			//numSubdomains.Add(new int[] { 128, 128, 128 });
			//numSubdomains.Add(new int[] { 256, 256, 256 });

			return (numElements, numSubdomains);
		}

		public (List<int[]> numElements, List<int[]> numSubdomains) GetParametricConfigConstSubdomainPerElementSize()
		{
			var numElements = new List<int[]>();
			var numSubdomains = new List<int[]>();

			numSubdomains.Add(new int[] { 3, 3, 3 });
			numSubdomains.Add(new int[] { 4, 4, 4 });
			numSubdomains.Add(new int[] { 5, 5, 5 });
			numSubdomains.Add(new int[] { 6, 6, 6 });
			numSubdomains.Add(new int[] { 7, 7, 7 });
			numSubdomains.Add(new int[] { 8, 8, 8 });
			numSubdomains.Add(new int[] { 9, 9, 9 });
			numSubdomains.Add(new int[] { 10, 10, 10 });
			numSubdomains.Add(new int[] { 11, 11, 11 });
			numSubdomains.Add(new int[] { 12, 12, 12 });
			numSubdomains.Add(new int[] { 13, 13, 13 });
			numSubdomains.Add(new int[] { 14, 14, 14 });
			numSubdomains.Add(new int[] { 15, 15, 15 });
			numSubdomains.Add(new int[] { 16, 16, 16 });
			numSubdomains.Add(new int[] { 17, 17, 17 });
			numSubdomains.Add(new int[] { 18, 18, 18 });
			numSubdomains.Add(new int[] { 19, 19, 19 });
			numSubdomains.Add(new int[] { 20, 20, 20 });

			numElements.Add(new int[] { 12, 12, 12 });
			numElements.Add(new int[] { 16, 16, 16 });
			numElements.Add(new int[] { 20, 20, 20 });
			numElements.Add(new int[] { 24, 24, 24 });
			numElements.Add(new int[] { 28, 28, 28 });
			numElements.Add(new int[] { 32, 32, 32 });
			numElements.Add(new int[] { 36, 36, 36 });
			numElements.Add(new int[] { 40, 40, 40 });
			numElements.Add(new int[] { 44, 44, 44 });
			numElements.Add(new int[] { 48, 48, 48 });
			numElements.Add(new int[] { 52, 52, 52 });
			numElements.Add(new int[] { 56, 56, 56 });
			numElements.Add(new int[] { 60, 60, 60 });
			numElements.Add(new int[] { 64, 64, 64 });
			numElements.Add(new int[] { 68, 68, 68 });
			numElements.Add(new int[] { 72, 72, 72 });
			numElements.Add(new int[] { 76, 76, 76 });
			numElements.Add(new int[] { 80, 80, 80 });

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

			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxX, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxX, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxX, StructuralDof.TranslationZ, 0.0);

			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinY, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinY, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinY, StructuralDof.TranslationZ, 0.0);

			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxY, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxY, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxY, StructuralDof.TranslationZ, 0.0);

			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinZ, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinZ, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MinZ, StructuralDof.TranslationZ, 0.0);

			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxZ, StructuralDof.TranslationX, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxZ, StructuralDof.TranslationY, 0.0);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.MaxZ, StructuralDof.TranslationZ, 0.0);

			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.Centroid, StructuralDof.TranslationX,
				CenterDisplacement[0]);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.Centroid, StructuralDof.TranslationY,
				CenterDisplacement[1]);
			builder.PrescribeDisplacement(UniformDdmModelBuilder3D.BoundaryRegion.Centroid, StructuralDof.TranslationZ,
				CenterDisplacement[2]);

			return builder;
		}
	}
}
