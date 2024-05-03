namespace MGroup.Solvers.DDM.DiscretizationExtensions
{
	using System.Collections.Concurrent;

	using MGroup.MSolve.DataStructures;

	public class TableExtension<TRow, TColumn, TValue> : Table<TRow, TColumn, TValue>
	{
		public bool TryRemove(TRow row, TColumn col)
		{
			bool containsRow = this.data.TryGetValue(row, out ConcurrentDictionary<TColumn, TValue> wholeRow);
			if (containsRow)
			{
				return wholeRow.TryRemove(col, out _);
			}
			else
			{
				return false;
			}
		}
	}
}
