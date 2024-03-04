using System;
using System.Collections.Generic;
using MGroup.LinearAlgebra.Distributed;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.Solvers.DDM.LinearSystem
{
	public class DistributedLinearSystem<TMatrix> : IGlobalLinearSystem
		where TMatrix : class, IMatrix
	{
		private readonly Func<IGlobalVector, DistributedOverlappingVector> checkCompatibleVector;
		private readonly Func<IGlobalMatrix, DistributedOverlappingMatrix<TMatrix>> checkCompatibleMatrix;

		public DistributedLinearSystem(Func<IGlobalVector, DistributedOverlappingVector> checkCompatibleVector,
			Func<IGlobalMatrix, DistributedOverlappingMatrix<TMatrix>> checkCompatibleMatrix)
		{
			this.checkCompatibleVector = checkCompatibleVector;
			this.checkCompatibleMatrix = checkCompatibleMatrix;
			Observers = new HashSet<ILinearSystemObserver>();
		}

		IGlobalMatrix IGlobalLinearSystem.Matrix
		{
			get => Matrix;
			set
			{
				DistributedOverlappingMatrix<TMatrix> globalMatrix = checkCompatibleMatrix(value);
				foreach (var observer in Observers)
				{
					observer.HandleMatrixWillBeSet();
				}
				Matrix = globalMatrix;
			}
		}

		public DistributedOverlappingMatrix<TMatrix> Matrix { get; set; }

		public HashSet<ILinearSystemObserver> Observers { get; }

		IGlobalVector IGlobalLinearSystem.RhsVector
		{
			get => RhsVector;
			set
			{
				DistributedOverlappingVector globalVector = checkCompatibleVector(value);
				RhsVector = globalVector;
			}
		}

		public DistributedOverlappingVector RhsVector { get; set; }

		IGlobalVector IGlobalLinearSystem.Solution => Solution;

		public DistributedOverlappingVector Solution { get; set; }
	}
}
