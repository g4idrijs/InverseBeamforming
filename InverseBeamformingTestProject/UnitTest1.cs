using InverseBeamforming;
using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace InverseBeamforming.Tests
{
	[TestClass()]
	public class UnitTest1
	{
		[TestMethod()]
		public void PositionsTest()
		{
			double[,] txpos = { { 1, 0, 0 }, { -1, 0, 0 } };
			double[,] rxpos = { { 10, 0, 0 } };

			Positions pos = new Positions(txpos, rxpos);

			Assert.AreEqual(txpos, pos.TxPositions);
			Assert.AreEqual(rxpos, pos.RxPositions);
		}

		[TestMethod()]
		public void getDistancesTest()
		{
			double[,] txpos = { { 1, 0, 0 }, { -1, 0, 0 } };
			double[,] rxpos = { { 0, 10, 0 } };

			Positions pos = new Positions(txpos, rxpos);

			double[] distances = { Math.Sqrt(101), Math.Sqrt(101) };

			var posDistances = pos.getDistances(0);

			for (int i=0; i<distances.Length; i++)
				Assert.AreEqual(distances[i], posDistances[i], 0.01);
		}
	}
}

