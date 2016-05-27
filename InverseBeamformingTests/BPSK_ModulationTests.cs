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
	public class BPSK_ModulationTests
	{
		[TestMethod()]
		public void BPSK_ModulationTest()
		{
			var bpsk = new BPSK_Modulation(10, 40, 1, 20, .5);

		}

		[TestMethod()]
		public void ModulateBitsTest()
		{
			var bpsk = new BPSK_Modulation(10, 100, 1, 20, .5);
			double[] waveform=bpsk.ModulateBits(10);
		}
	}
}