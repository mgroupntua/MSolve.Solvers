using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public abstract class SubmatrixExtractorCsrBase
	{
		protected CsrMatrix originalMatrix;

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
		/// Maps indices of the values arrays between the original matrix A and its submatrix A01
		/// A10.Values[i] = A.Values[map10[i]]
		/// </summary>
		protected int[] map10;

		/// <summary>
		/// Maps indices of the values arrays between the original matrix A and its submatrix A11
		/// A11.Values[i] = A.Values[map11[i]]
		/// </summary>
		protected int[] map11;

		public SubmatrixExtractorCsrBase()
		{
		}

		public CsrMatrix Submatrix01 { get; protected set; }

		public CsrMatrix Submatrix10 { get; protected set; }

		public CscMatrix Submatrix11 { get; protected set; }

		public virtual void Clear()
		{
			originalMatrix = null;
			Submatrix01 = null;
			Submatrix10 = null;
			Submatrix11 = null;
			map00 = null;
			map01 = null;
			map10 = null;
			map11 = null;
		}
	}
}
