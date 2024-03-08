using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.FetiDP.StiffnessMatrices
{
	public interface IFetiDPSubdomainMatrixManager
	{
		bool IsEmpty { get; }

		IMatrix SchurComplementOfRemainderDofs { get; }

		void CalcSchurComplementOfRemainderDofs();

		void ClearSubMatrices();

		/// <summary>
		/// Divide subdomain matrix between internal and boundary-remainder dofs.
		/// </summary>
		void ExtractKiiKbbKib();

		void ExtractKrrKccKrc();

		void HandleDofsWereModified();

		void InvertKii(bool diagonalOnly);

		void InvertKrr();

		Vector MultiplyInverseKiiTimes(Vector vector, bool diagonalOnly);

		Vector MultiplyInverseKrrTimes(Vector vector);

		Vector MultiplyKbbTimes(Vector vector);

		Vector MultiplyKbiTimes(Vector vector);

		Vector MultiplyKccTimes(Vector vector);

		Vector MultiplyKcrTimes(Vector vector);

		Vector MultiplyKibTimes(Vector vector);

		Vector MultiplyKrcTimes(Vector vector);

		/// <summary>
		/// S[s] * x = (Kcc[s] - Kcr[s] * inv(Krr[s]) * Krc[s]) * x, where s is a subdomain and x is the <paramref name="input"/>.
		/// </summary>
		/// <param name="input">The displacements that correspond to corner dofs of this subdomain.</param>
		/// <param name="output">The forces that correspond to corner dofs of this subdomain.</param>
		public void MultiplySchurComplementImplicitly(Vector input, Vector output)
		{
			output.CopyFrom(MultiplyKccTimes(input));
			Vector temp = MultiplyKrcTimes(input);
			temp = MultiplyInverseKrrTimes(temp);
			temp = MultiplyKcrTimes(temp);
			output.SubtractIntoThis(temp);
		}

		void ReorderInternalDofs();

		void ReorderRemainderDofs();

		Matrix CalcInvKrrTimesKrc();
	}
}
