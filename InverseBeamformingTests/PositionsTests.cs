using Microsoft.VisualStudio.TestTools.UnitTesting;
using InverseBeamforming;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming.Tests
{

	[TestClass()]
	public class PositionsTests
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

			for (int i = 0; i < distances.Length; i++)
				Assert.AreEqual(distances[i], posDistances[i], 0.01);
		}

		[TestMethod()]
		public void getAllDistancesTest()
		{
			double[,] txpos = { { 1, 0, 0 }, { -1, 0, 0 } };
			double[,] rxpos = { { 0, 10, 0 }, { 10, 100, 0 } };

			Positions pos = new Positions(txpos, rxpos);

			double[,] distances = { { Math.Sqrt(101), Math.Sqrt(10081) }, { Math.Sqrt(101) , Math.Sqrt(10121) } };

			var posDistances = pos.getAllDistances();

			for (int i = 0; i < 2; i++)
				for(int k=0; k<2; k++)
					Assert.AreEqual(distances[i,k], posDistances[i,k], 0.01);
		}
	}
}
