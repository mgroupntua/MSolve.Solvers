using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

//TODOMPI: Standardize the mappings. I think 100% of the cases, the following case holds: 
//		There are 2 sets of entries, the full entries and a subset. There is an array that takes indices from the subset and
//		its values are in the full entries. E.g. localToGlobal[i] = j, where i = 0 ... numLocal-1 and j = 0 ... numGlobal-1.
//		In some cases there is also scaling represented by different arrays. Usually clients need to multiply vectors 
//		with this mapping matrix or its transpose. In some cases, clients need to access the intrnal mapping array (e.g. during 
//		FEM assembly). Both of these functionalities should be provided. I should also be consistent in using them. Some clients
//		may expose one or the other concrete class or even only the internal array. Perhaps it would be better to always expose  
//		this interface (e.g. DofSeparator components should expose IMappingMatrix. MatrixManager components should take the 
//		internal array exposed from the IMappingMatrix).
namespace MGroup.Solvers.DDM.Mappings
{
	public interface IMappingMatrix
	{
		int NumColumns { get; }

		int NumRows { get; }

		Matrix CopyToFullMatrix();

		Vector Multiply(Vector vector, bool transpose);
	}
}
