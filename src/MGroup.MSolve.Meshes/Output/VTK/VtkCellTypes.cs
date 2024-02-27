using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization;

//TODOMesh: allow reordering e.g. for Hexa8
namespace MGroup.MSolve.Meshes.Output.VTK
{
	public static class VtkCellTypes
	{
		private static readonly IReadOnlyDictionary<CellType, int> cellTypeCodes = 
			new Dictionary<CellType, int>
			{
				{ CellType.Line2, 3 },

				//{ CellType.Polygon, 7 },

											// 3 ---- 2
											// |      |
											// |      |
				{ CellType.Quad4, 9 },      // 0 ---- 1

											// 3 -- 6 -- 2
											// |         |
											// 7         5
											// |         |
				{ CellType.Quad8, 23 },     // 0 -- 4 -- 1

											// 3 -- 6 -- 2
											// |    |    |
											// 7 -- 8 -- 5
											// |    |    |
				{ CellType.Quad9, 28 },     // 0 -- 4 -- 1

											//    2
											//   /  \
											//  /    \
				{ CellType.Tri3, 5 },       // 0 ---  1

											//     2
											//    /  \
											//   5    4
											//  /       \
				{ CellType.Tri6, 22 },       // 0 -- 3 -- 1

				// 3----------2
				// |\         |\
				// | \        | \
				// |  \       |  \
				// |   7------+---6
				// |   |      |   |
				// 0---+------1   |
				//  \  |       \  |
				//   \ |        \ |
				//    \|         \|
				//     4----------5
				{ CellType.Hexa8, 12 },

				{CellType.Tet4, 10 }
			};

		public static CellType GetMSolveCellType(int vtkCellCode) => cellTypeCodes.First(pair => pair.Value == vtkCellCode).Key;

		public static int GetVtkCellCode(CellType msolveCellType) => cellTypeCodes[msolveCellType];
	}
}
