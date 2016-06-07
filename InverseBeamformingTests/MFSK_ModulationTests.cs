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
	public class MFSK_ModulationTests
	{
		private string _testFileDumpDirec = @"C:\Users\vikin_000\OneDrive\AFIT\Thesis\Code\C_Sharp\InvBeamLib\InverseBeamforming\TestFileDump\FSK" + Path.DirectorySeparatorChar;

		[TestMethod()]
		public void MFSK_ModulationTest()
		{
			var mfsk = new MFSK_Modulation(10, 40, 1, 20, .5, null, 1000, 4, new double[] { 1000, 2000, 3000, 4000 });
		}

		[TestMethod()]
		public void ModulateBitsTest()
		{
			var sampleRate = 1000;
			var mfsk = new MFSK_Modulation(10, 40, sampleRate, 100, .5, null, 1000, 4, new double[] { 10, 20, 30, 40 });
			double[] waveform = mfsk.ModulateBits();
			writeToCSV(waveform, "FSK_ModulateBitsTest");
		}


		[TestMethod()]
		public void DemodulateWaveformTest()
		{
			var sampleRate = 1000;
			var mfsk = new MFSK_Modulation(10, 40, sampleRate, 100, .5, null, 1000, 4, new double[] { 10, 20, 30, 40 });

			byte[] inbits = mfsk.GenerateRandomBits();
			double[] waveform = mfsk.ModulateBits(inbits);
			byte[] outbits = mfsk.CorrelationReceiver(waveform);
			writeToCSV(waveform, "FSK_Modulation",sampleRate);
			writeToCSV(inbits, "inbits");
			writeToCSV(outbits, "outbits");
			CollectionAssert.AreEqual(inbits, outbits);
		}

		
		private void writeToCSV<T>(T[] mat, string filename,double sampleRate=0)
		{
			filename = _testFileDumpDirec + filename + ".csv";
			var csv = new StringBuilder();

			if (sampleRate!=0)
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