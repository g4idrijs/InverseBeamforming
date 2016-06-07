using Microsoft.VisualStudio.TestTools.UnitTesting;
using InverseBeamforming;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace InverseBeamforming.Tests
{
	[TestClass()]
	public class MPSK_ModulationTests
	{
		private string _testFileDumpDirec = @"C:\Users\vikin_000\OneDrive\AFIT\Thesis\Code\C_Sharp\InvBeamLib\InverseBeamforming\TestFileDump\MPSK" + Path.DirectorySeparatorChar;

		[TestMethod()]
		public void MPSK_ModulationTest()
		{
			var mpsk = new MPSK_Modulation(10, 40, 1, 20, .5, null, 1000, 4, new double[] { 0, Math.PI/4, Math.PI / 2, 3*Math.PI / 4 });
		}

		[TestMethod()]
		public void ModulateBitsTest()
		{
			var sampleRate = 1000;
			var mpsk = new MPSK_Modulation(100, sampleRate, 40, 100, .5, null, 1000, 4, new double[] { 0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4 });
			double[] waveform = mpsk.ModulateBits();
			writeToCSV(waveform, "MPSK_ModulateBitsTest");
		}

		[TestMethod()]
		public void DemodulateWaveformTest()
		{
			var sampleRate = 1000;
			var mpsk = new MPSK_Modulation(100, sampleRate, 40, 100, .5, null, 1000, 4, new double[] { 0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4 });
			byte[] inbits = mpsk.GenerateRandomBits();
			double[] waveform = mpsk.ModulateBits(inbits);
			byte[] outbits = mpsk.CorrelationReceiver(waveform);
			writeToCSV(waveform, "MPSK_Modulation", sampleRate);
			writeToCSV(inbits, "inbits");
			writeToCSV(outbits, "outbits");
			CollectionAssert.AreEqual(inbits, outbits);
		}

		private void writeToCSV<T>(T[] mat, string filename, double sampleRate = 0)
		{
			filename = _testFileDumpDirec + filename + ".csv";
			var csv = new StringBuilder();

			if (sampleRate != 0)
				csv.AppendLine(sampleRate.ToString());

			for (int i = 0; i < mat.Length - 1; i++)
			{
				csv.AppendLine(mat[i].ToString() + ",");
			}
			csv.AppendLine(mat[mat.Length - 1].ToString());

			File.WriteAllText(filename, csv.ToString());
		}
	}
}