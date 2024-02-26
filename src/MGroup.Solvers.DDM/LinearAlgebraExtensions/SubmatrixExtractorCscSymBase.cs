using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public abstract class SubmatrixExtractorCscSymBase
	{
		protected SymmetricCscMatrix originalMatrix;

		/// <summary>
		/// Maps indices of the values arrays between the original matrix A and its submatrix A00
		/// A00.Values[i] = A.Values[map00[i]]
		/// </summary>
		protected int[] map00;

		/// <summary>
		/// Maps indices of the values arrays between the original matrix A and its submatrix A01
		/// A01.Values[i] = A.Values[map01[i]]
		/// </summary>
		protected int[] map01;

		/// <summary>
		/// Maps indices of the values arrays between the original matrix A and its submatrix A11
		/// A11.Values[i] = A.Values[map11[i]]
		/// </summary>
		protected int[] map11;

		public SubmatrixExtractorCscSymBase()
		{
		}

		public CsrMatrix Submatrix01 { get; protected set; }

		public SymmetricCscMatrix Submatrix11 { get; protected set; }

		public virtual void Clear()
		{
			originalMatrix = null;
			Submatrix01 = null;
			Submatrix11 = null;
			map00 = null;
			map01 = null;
			map11 = null;
		}

		/// <summary>
		/// The sparsity pattern in symmetric CSC format of a submatrix K00, with <paramref name="indicesGroup0"/> specifying 
		/// which rows and columns to keep.
		/// </summary>
		public (int[] rowIndices, int[] colOffsers) ExtractSparsityPattern(SymmetricCscMatrix originalMatrix, int[] indicesGroup0)
		{
			// Store which entries of each column are nonzero
			var columnsA00 = new SortedSet<int>[indicesGroup0.Length];
			for (int j = 0; j < indicesGroup0.Length; ++j)
			{
				columnsA00[j] = new SortedSet<int>();
			}

			// Original matrix indices to submatrix indices. Indices not belonging to group 0 will be marked as -1.
			var originalToSubIndices = new int[originalMatrix.NumRows];
			Fill(originalToSubIndices, -1);
			for (int i0 = 0; i0 < indicesGroup0.Length; ++i0)
			{
				originalToSubIndices[indicesGroup0[i0]] = -i0 - 1;
			}

			// Iterate the non zero values array of the original sparse matrix
			for (int j = 0; j < originalMatrix.NumColumns; ++j)
			{
				int start = originalMatrix.RawColOffsets[j];
				int end = originalMatrix.RawColOffsets[j + 1];
				int subJ = originalToSubIndices[j];
				if (subJ >= 0) // j belongs to group 0 indices
				{
					for (int t = start; t < end; ++t)
					{
						int i = originalMatrix.RawRowIndices[t];
						int subI = originalToSubIndices[i];
						if (subI >= 0) // (i, j) belongs to submatrix A00
						{
							columnsA00[subJ].Add(subI);
						}
					}
				}
			}

			return BuildCscArrays(indicesGroup0.Length, columnsA00);
		}

		public static (int[] rowIndices, int[] colOffsers) BuildCscArrays(int order, SortedSet<int>[] nonzeroRowsOfEachCol)
		{
			// Create CSC arrays from the dictionary
			int[] colOffsets = new int[order + 1];
			int nnz = 0;
			for (int j = 0; j < order; ++j)
			{
				colOffsets[j] = nnz;
				nnz += nonzeroRowsOfEachCol[j].Count;
			}
			colOffsets[order] = nnz; //The last CSC entry is nnz.

			int[] rowIndices = new int[nnz];
			double[] values = new double[nnz];
			int counter = 0;
			for (int j = 0; j < order; ++j)
			{
				foreach (var rowIdx in nonzeroRowsOfEachCol[j])
				{
					rowIndices[counter] = rowIdx;
					++counter;
				}
			}

			return (rowIndices, colOffsets);
		}

		/// <summary>
		/// <paramref name="submatrixValues"/>[i] = <paramref name="originalValues"/>[<paramref name="map"/>[i]]
		/// </summary>
		/// <param name="originalValues"></param>
		/// <param name="submatrixValues"></param>
		/// <param name="map">All entries correspond to indices into <paramref name="originalValues"/>.</param>
		public static void CopyValuesArray(double[] originalValues, double[] submatrixValues, int[] map)
		{
			for (int i = 0; i < submatrixValues.Length; ++i)
			{
				submatrixValues[i] = originalValues[map[i]];
			}
		}

		/// <summary>
		/// Same as <see cref="CopyValuesArray(double[], double[], int[])"/>, but takes into account that some entries of
		/// <paramref name="submatrixValues"/> are zeros, not originally stored in <paramref name="originalValues"/>.
		/// If <paramref name="map"/>[i] negative then <paramref name="submatrixValues"/>[i] = 0, else
		/// <paramref name="submatrixValues"/>[i] = <paramref name="originalValues"/>[<paramref name="map"/>[i]]
		/// </summary>
		/// <param name="originalValues"></param>
		/// <param name="submatrixValues"></param>
		/// <param name="map">All entries correspond to indices into <paramref name="originalValues"/>.</param>
		public static void CopyValuesArrayAndZeros(double[] originalValues, double[] submatrixValues, int[] map)
		{
			for (int i = 0; i < submatrixValues.Length; ++i)
			{
				int indexOriginal = map[i];
				if (indexOriginal >= 0)
				{
					submatrixValues[i] = originalValues[indexOriginal];
				}
				else
				{
					submatrixValues[i] = 0.0;
				}
			}
		}

		//TODO: There must be a faster way to do this
		public static void Fill<T>(T[] array, T value)
		{
			for (int i = 0; i < array.Length; ++i)
			{
				array[i] = value;
			}
		}

		/// <summary>
		/// Creates an array that maps the indices (rows/columns) of the original matrix to the indices of its 2x2 submatrices. 
		/// Group 0 indices i are stored as: -i-1. Group 1 indices are stored as they are.
		/// </summary>
		/// <param name="originalOrder"></param>
		/// <param name="indicesGroup0"></param>
		/// <param name="indicesGroup1"></param>
		public static int[] MapOriginalToSubmatrixIndices(int originalOrder, int[] indicesGroup0, int[] indicesGroup1)
		{
			var originalToSubIndices = new int[originalOrder];
			for (int i0 = 0; i0 < indicesGroup0.Length; ++i0)
			{
				originalToSubIndices[indicesGroup0[i0]] = -i0 - 1;
			}

			for (int i1 = 0; i1 < indicesGroup1.Length; ++i1)
			{
				originalToSubIndices[indicesGroup1[i1]] = i1;
			}

			return originalToSubIndices;
		}
	}
}
