namespace MGroup.Solvers.DDM.LinearAlgebraExtensions.Matrices
{
	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public static class MatrixExtensions
	{
		/// <summary>
		/// Returns MGroup.LinearAlgebra.Vectors.Vector with the entries of the matrix's main diagonal.
		/// </summary>
		/// <param name="matrix">The matrix whose diagonal will be copied. It must be square.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public static Vector GetDiagonal(this IMatrixView matrix)
		{
			return Vector.CreateFromArray(matrix.GetDiagonalAsArray());
		}

		/// <summary>
		/// Returns an array with the entries of the matrix's main diagonal.
		/// </summary>
		/// <param name="matrix">The matrix whose diagonal will be copied. It must be square.</param>
		/// <exception cref="NonMatchingDimensionsException">Thrown if the matrix is not square.</exception>
		public static double[] GetDiagonalAsArray(this IMatrixView matrix)
		{
			Preconditions.CheckSquare(matrix);
			double[] array = new double[matrix.NumRows];
			for (int i = 0; i < matrix.NumRows; i++)
			{
				array[i] = matrix[i, i];
			}

			return array;
		}
	}
}
