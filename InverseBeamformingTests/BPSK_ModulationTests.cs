using Microsoft.VisualStudio.TestTools.UnitTesting;
using InverseBeamforming;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace InverseBeamforming.Tests
{
	[TestClass()]
	public class BPSK_ModulationTests
	{
		[TestMethod()]
		public void BPSK_ModulationTest()
		{
			var bpsk = new BPSK_Modulation(10, 40, 1, 20, .5, null);

		}

		[TestMethod()]
		public void ModulateBitsTest()
		{
			var bpsk = new BPSK_Modulation(10, 100, 1, 20, .5, null);
			double[] waveform = bpsk.ModulateBits(10);
		}

		[TestMethod()]
		public void FIR_FilterTest()
		{
			double[] coefs = { 0.0004918999, 0.0003199292, -0.0000000000, -0.0003633909, -0.0006312609, -0.0006698830, -0.0004061822, 0.0001195890, 0.0007335821, 0.0011668717, 0.0011566730, 0.0005795375, -0.0004436998, -0.0015362768, -0.0021818714, -0.0019451886, -0.0007111603, 0.0011751752, 0.0029643583, 0.0037731036, 0.0029846847, 0.0006097011, -0.0025792523, -0.0052531992, -0.0060437038, -0.0041851771, 0.0000000000, 0.0050347924, 0.0087507230, 0.0091632713, 0.0054278953, -0.0015509929, -0.0092006763, -0.0141371985, -0.0135491320, -0.0065803763, 0.0049027331, 0.0166053238, 0.0232167449, 0.0205322359, 0.0075140558, -0.0125643771, -0.0324946486, -0.0431063759, -0.0362982920, -0.0081217790, 0.0392103634, 0.0972196347, 0.1531611401, 0.1935823237, 0.2083155053, 0.1935823237, 0.1531611401, 0.0972196347, 0.0392103634, -0.0081217790, -0.0362982920, -0.0431063759, -0.0324946486, -0.0125643771, 0.0075140558, 0.0205322359, 0.0232167449, 0.0166053238, 0.0049027331, -0.0065803763, -0.0135491320, -0.0141371985, -0.0092006763, -0.0015509929, 0.0054278953, 0.0091632713, 0.0087507230, 0.0050347924, 0.0000000000, -0.0041851771, -0.0060437038, -0.0052531992, -0.0025792523, 0.0006097011, 0.0029846847, 0.0037731036, 0.0029643583, 0.0011751752, -0.0007111603, -0.0019451886, -0.0021818714, -0.0015362768, -0.0004436998, 0.0005795375, 0.0011566730, 0.0011668717, 0.0007335821, 0.0001195890, -0.0004061822, -0.0006698830, -0.0006312609, -0.0003633909, -0.0000000000, 0.0003199292, 0.0004918999 };
			var bpsk = new BPSK_Modulation(10, 100, 1, 20, .5, coefs);
			double[] signal = new double[1024];
			double f1 = 3000, f2 = 8000, fs = 48000;
			double pi2time;

			for (int i = 0; i < 1024; i++)
			{
				pi2time = 2 * Math.PI * (double)i / fs;
				signal[i] = Math.Sin(pi2time * f1) + Math.Sin(pi2time * f2);
			}

			signal = bpsk.FIR_Filter(signal);
		}

		[TestMethod()]
		public void DemodulateWaveformTest()
		{
			var bpsk = new BPSK_Modulation(10, 100, 1, 20, .5, null);

			byte[] inbits = bpsk.GenerateRandomBits(1000);
			double[] waveform = bpsk.ModulateBits(inbits);
			byte[] outbits = bpsk.DemodulateWaveform(waveform);

			CollectionAssert.AreEqual(inbits, outbits);
		}

		[TestMethod()]
		public void AdditiveWhiteGaussianNoiseTest()
		{
			var bpsk = new BPSK_Modulation(10, 100, 1, 20, .5, null);
			double[] waveform = bpsk.AdditiveWhiteGaussianNoise(10000,1);
		}

		[TestMethod()]
		public void numberDifferentBitsTest()
		{
			var bpsk = new BPSK_Modulation(10, 100, 1, 20, .5, null);
			byte[] inbits = { 1, 0, 1, 1, 0, 1 };
			byte[] outbits = { 1, 1, 0, 0, 0, 1 };

			Assert.AreEqual(3, bpsk.numberDifferentBits(inbits, outbits));
			CollectionAssert.AreEqual(new int[] { 1, 2 }, bpsk.numberDifferentBitsEachType(inbits, outbits));
		}

		[TestMethod()]
		public void Mod_and_DemodManyBitsTest()
		{
			var bpsk = new BPSK_Modulation(10, 100, 1, 20, 1, null);

			int numTests = 10;
			byte[] inbits = bpsk.GenerateRandomBits(1000);
			double[] waveform = bpsk.ModulateBits(inbits);
			double[] noisePower = new double[numTests];
			byte[] outbits = bpsk.DemodulateWaveform(waveform);

			int[,] numWrong = new int[numTests, 2];
			int[] temp;
			for(int i=0; i< numTests; i++)
			{
				inbits = bpsk.GenerateRandomBits(1000);
				waveform = bpsk.ModulateBits(inbits);
				noisePower[i] = bpsk.AdditiveWhiteGaussianNoise(ref waveform, 1);
				outbits = bpsk.DemodulateWaveform(waveform);
				temp = bpsk.numberDifferentBitsEachType(inbits, outbits);
				numWrong[i, 0] = temp[0];
				numWrong[i, 1] = temp[1];
			}
		}

		[TestMethod()]
		public void RunSimulation_Test()
		{
			var bpsk = new BPSK_Modulation(101680, 1220160, 1, 7680, .5, null);
			double ber=bpsk.RunSimulationOneNoisePower(500, 1000, 2000);
			Debug.Print(ber.ToString());
		}
	}
}