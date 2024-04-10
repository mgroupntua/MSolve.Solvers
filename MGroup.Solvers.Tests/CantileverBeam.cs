// Geometry:
//             | 
//             |
//            \|/ load
// |           V 
// |------------
// |    length

// Cross section:
//     |
//     | 
//    \|/ load
//     .
//   -----
//  |     |
//  |     |
//  |     | <--- height
//  |     |
//  |     |
//   -----
//     ^
//     |     
//   width

namespace MGroup.Solvers.Tests
{
	using MGroup.Constitutive.Structural;
	using MGroup.MSolve.Solution.LinearSystem;
	using MGroup.MSolve.Solution.AlgebraicModel;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.Constitutive.Structural.Planar;
	using MGroup.MSolve.Discretization.Meshes.Generation.Custom;
	using MGroup.MSolve.Discretization.Meshes.Generation;
	using MGroup.FEM.Structural.Continuum;
	using MGroup.Constitutive.Structural.BoundaryConditions;
	using MGroup.Constitutive.Structural.Transient;

	public class CantileverBeam
	{
		private const double eulerBeamDimensionsRatio = 10.0;
		private const int subdomainID = 0;

		private readonly double length;
		private readonly double height;
		private readonly double width;
		private readonly double endPointLoad;
		private readonly double youngModulus;
		private readonly INode[] endNodes;

		private CantileverBeam(double length, double height, double width, double endPointLoad,
			double youngModulus, Model model, INode[] endNodes)
		{
			this.length = length;
			this.height = height;
			this.width = width;
			this.endPointLoad = endPointLoad;
			this.youngModulus = youngModulus;
			this.Model = model;
			this.endNodes = endNodes;
		}

		public Model Model { get; }

		public double CalculateEndDeflectionWithEulerBeamTheory()
		{
			if ((length < eulerBeamDimensionsRatio * height) || (length < eulerBeamDimensionsRatio * width))
			{
				throw new ArgumentException("In order for Euler beam theory to hold, the beam's length must be at least"
					+ $" {eulerBeamDimensionsRatio} times greater than each of the other two dimensions.");
			}
			double momentOfInertia = 1.0 / 12.0 * width * Math.Pow(height, 3.0);
			return endPointLoad * Math.Pow(length, 3.0) / (3.0 * youngModulus * momentOfInertia);
		}

		public double CalculateAverageEndDeflectionFromSolution(IGlobalVector solution, IVectorValueExtractor resultsExtractor)
		{
			//DofTable subdomainDofs = solver.LinearSystem.DofOrdering.SubdomainDofOrderings[subdomainID].FreeDofs;
			double endDeflectionSum = 0.0;
			int dofsCount = 0;
			foreach (INode node in endNodes)
			{
				try
				{
					double displacement = resultsExtractor.ExtractSingleValue(solution, node, StructuralDof.TranslationY);
					++dofsCount;
					endDeflectionSum += displacement;
				}
				catch (KeyNotFoundException)
				{
					// Ignore constrained dofs
				}
			}
			return endDeflectionSum / dofsCount;
		}


		public class Builder
		{
			public double Length { get; set; } = 10.00;
			public double Height { get; set; } = 0.50;
			public double Width { get; set; } = 0.25;
			public double EndPointLoad { get; set; } = 20.0E3;
			public double YoungModulus { get; set; } = 2.1E7;
			public double PoissonRatio { get; set; } = 0.3;

			public CantileverBeam BuildWithBeamElements(int numElements)
			{
				throw new NotImplementedException();
			}

			public CantileverBeam BuildWithQuad4Elements(int numElementsAlongLength, int numElementsAlongHeight)
			{
				// Material and section properties
				double thickness = Width;
				var material = new ElasticMaterial2D(YoungModulus, PoissonRatio, StressState2D.PlaneStress);
				var dynamicProperties = new TransientAnalysisProperties(1.0, 0.0, 0.0);

				// Model with 1 subdomain
				var model = new Model();
				model.SubdomainsDictionary.Add(subdomainID, new Subdomain(subdomainID));

				// Generate mesh
				var meshGenerator = new UniformMeshGenerator2D<Node>(0.0, 0.0, Length, Height,
					numElementsAlongLength, numElementsAlongHeight);
				(IReadOnlyList<Node> vertices, IReadOnlyList<CellConnectivity<Node>> cells) =
					meshGenerator.CreateMesh((id, x, y, z) => new Node(id: id, x: x, y: y, z: z));

				// Add nodes to the model
				for (int n = 0; n < vertices.Count; ++n) model.NodesDictionary.Add(n, vertices[n]);

				// Add Quad4 elements to the model
				var factory = new ContinuumElement2DFactory(thickness, material, dynamicProperties);
				for (int e = 0; e < cells.Count; ++e)
				{
					ContinuumElement2D element = factory.CreateElement(cells[e].CellType, cells[e].Vertices);
					element.ID = e;
					model.ElementsDictionary.Add(e, element);
					model.SubdomainsDictionary[subdomainID].Elements.Add(element);
				}

				// Clamp boundary condition at one end
				double tol = 1E-10; //TODO: this should be chosen w.r.t. the element size along X
				var dirichletBCs = new List<INodalDisplacementBoundaryCondition>();
				foreach (INode node in model.NodesDictionary.Values.Where(node => Math.Abs(node.X) <= tol))
				{
					dirichletBCs.Add(new NodalDisplacement(node, StructuralDof.TranslationX, 0));
					dirichletBCs.Add(new NodalDisplacement(node, StructuralDof.TranslationY, 0));
				}

				// Apply concentrated load at the other end
				INode[] loadedNodes = model.NodesDictionary.Values.Where(node => Math.Abs(node.X - Length) <= tol).ToArray();
				double loadPerNode = EndPointLoad / loadedNodes.Length;
				var neumannBCs = new List<INodalLoadBoundaryCondition>();
				foreach (INode node in loadedNodes)
				{
					neumannBCs.Add(new NodalLoad(node, StructuralDof.TranslationY, loadPerNode));
				}

				model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(dirichletBCs, neumannBCs));

				return new CantileverBeam(Length, Height, Width, EndPointLoad, YoungModulus, model, loadedNodes);
			}

			public CantileverBeam BuildWithHexa8Elements(int numElementsAlongLength, int numElementsAlongHeight,
				int numElementsAlongWidth)
			{
				throw new NotImplementedException();
			}
		}
	}
}
