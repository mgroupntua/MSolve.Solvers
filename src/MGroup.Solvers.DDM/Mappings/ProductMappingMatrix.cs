using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.Mappings
{
	public class ProductMappingMatrix : IMappingMatrix
	{
		private readonly IList<IMappingMatrix> matricesLeftToRight;
		public ProductMappingMatrix(IList<IMappingMatrix> matricesLeftToRight)
		{
			this.matricesLeftToRight = matricesLeftToRight;
			NumRows = matricesLeftToRight[0].NumRows;
			NumColumns = matricesLeftToRight.Last().NumColumns;
		}

		public int NumColumns { get; }

		public int NumRows { get; }

		public Matrix CopyToFullMatrix()
		{
			Matrix result = matricesLeftToRight[0].CopyToFullMatrix();
			for (int i = 1; i < matricesLeftToRight.Count; ++i)
			{
				result *= matricesLeftToRight[i].CopyToFullMatrix();
			}
			return result;
		}

		public Vector Multiply(Vector vector, bool transpose)
		{
			Vector temp = vector;
			if (transpose)
			{
				for (int i = 0; i < matricesLeftToRight.Count; i++)
				{
					temp = matricesLeftToRight[i].Multiply(temp, true);
				}
				return temp;
			}
			else
			{
				for (int i = matricesLeftToRight.Count - 1; i >= 0; i--)
				{
					temp = matricesLeftToRight[i].Multiply(temp, false);
				}
				return temp;
			}
		}
	}
}
