using Microsoft.VisualStudio.TestTools.UnitTesting;
using InverseBeamforming;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MPSK_Modulation=InverseBeamforming.Modulations.MPSK_Modulation;

namespace InverseBeamforming.Tests
{
	[TestClass()]
	public class MPSK_ModulationTests
	{
		private string _testFileDumpDirec = @"C:\Users\vikin\OneDrive\AFIT\Thesis\Code\C_Sharp\InvBeamLib\InverseBeamforming\TestFileDump\MPSK" + Path.DirectorySeparatorChar;

		[TestMethod()]
		public void MPSK_ModulationTest()
		{
			var mpsk = new MPSK_Modulation(10, 40, 1, 20, .5, null, 1000, 4, new double[] { 0, Math.PI/4, Math.PI / 2, 3*Math.PI / 4 }, new byte[] { 0, 1, 2, 3 });
		}

		[TestMethod()]
		public void ModulateBitsTest()
		{
			var sampleRate = 1000;
			var mpsk = new MPSK_Modulation(100, sampleRate, 40, 100, .5, null, 1000, 4, new double[] { 0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4 }, new byte[] { 0, 1, 2, 3 });
			double[] waveform = mpsk.ModulateBits();
			writeToCSV(waveform, "MPSK_ModulateBitsTest");
		}

		[TestMethod()]
		public void ModulateBits_GenerateWaveformFor_qa_ModulateBits_char_float_Test()
		{
			var seed = 40;
			var sampleRate = 100;
			var carrierFrequency = 10;
			var samplesPerSymbol = 10;
			var signalPower = 1;


			var mpsk = new MPSK_Modulation(carrierFrequency, sampleRate, seed, samplesPerSymbol, signalPower, null, 2, 4, new double[] { 0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4 }, new byte[] { 0, 1, 2, 3 });
			double[] waveform = mpsk.ModulateBits(new byte[] { 0, 1 });
			writeToCSV(waveform, "MPSK_ModulateBitsTest");
		}

		[TestMethod()]
		public void DemodulateWaveformTest()
		{
			var sampleRate = 1000;
			var mpsk = new MPSK_Modulation(100, sampleRate, 40, 100, .5, null, 1000, 4, new double[] { 0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4 }, new byte[] { 0, 1, 2, 3 });
			byte[] inbits = mpsk.GenerateRandomBits();
			double[] waveform = mpsk.ModulateBits(inbits);
			byte[] outbits = mpsk.CorrelationReceiver(waveform);
			writeToCSV(waveform, "MPSK_Modulation", sampleRate);
			writeToCSV(inbits, "inbits");
			writeToCSV(outbits, "outbits");
			CollectionAssert.AreEqual(inbits, outbits);
		}

		[TestMethod()]
		public void RunSimpleSimulationObservableTest()
		{
			var seed = 40;
			var sampleRate = 1220160;
			var carrierFrequency = 101680;
			var samplesPerSymbol = 7680;
			var signalPower = 1;
			var mpsk = new MPSK_Modulation(carrierFrequency, sampleRate, seed, samplesPerSymbol, signalPower, null, 2, 4, new double[] { 0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4 }, new byte[] { 0, 1, 2, 3 });
			

			string filename = "RunSimpleSimulationObservableTest";
			filename = _testFileDumpDirec + filename + ".csv";
			var obs = mpsk.Subscribe(new MPSK_Modulation.IntermediateSimResultsReporter(filename));
			mpsk.RunSimpleSimulationObservable(500, 1700);
			obs.Dispose();
		}
		
		[TestMethod()]
		public void RunFIRFilterSimulationObservableTest()
		{
			double[] coefs = { 0.0004918999, 0.0003199292, -0.0000000000, -0.0003633909, -0.0006312609, -0.0006698830, -0.0004061822, 0.0001195890, 0.0007335821, 0.0011668717, 0.0011566730, 0.0005795375, -0.0004436998, -0.0015362768, -0.0021818714, -0.0019451886, -0.0007111603, 0.0011751752, 0.0029643583, 0.0037731036, 0.0029846847, 0.0006097011, -0.0025792523, -0.0052531992, -0.0060437038, -0.0041851771, 0.0000000000, 0.0050347924, 0.0087507230, 0.0091632713, 0.0054278953, -0.0015509929, -0.0092006763, -0.0141371985, -0.0135491320, -0.0065803763, 0.0049027331, 0.0166053238, 0.0232167449, 0.0205322359, 0.0075140558, -0.0125643771, -0.0324946486, -0.0431063759, -0.0362982920, -0.0081217790, 0.0392103634, 0.0972196347, 0.1531611401, 0.1935823237, 0.2083155053, 0.1935823237, 0.1531611401, 0.0972196347, 0.0392103634, -0.0081217790, -0.0362982920, -0.0431063759, -0.0324946486, -0.0125643771, 0.0075140558, 0.0205322359, 0.0232167449, 0.0166053238, 0.0049027331, -0.0065803763, -0.0135491320, -0.0141371985, -0.0092006763, -0.0015509929, 0.0054278953, 0.0091632713, 0.0087507230, 0.0050347924, 0.0000000000, -0.0041851771, -0.0060437038, -0.0052531992, -0.0025792523, 0.0006097011, 0.0029846847, 0.0037731036, 0.0029643583, 0.0011751752, -0.0007111603, -0.0019451886, -0.0021818714, -0.0015362768, -0.0004436998, 0.0005795375, 0.0011566730, 0.0011668717, 0.0007335821, 0.0001195890, -0.0004061822, -0.0006698830, -0.0006312609, -0.0003633909, -0.0000000000, 0.0003199292, 0.0004918999 };
			var seed = 40;
			var sampleRate = 1220160;
			var carrierFrequency = 101680;
			var samplesPerSymbol = 7680;
			var signalPower = 1;
			var mpsk = new MPSK_Modulation(carrierFrequency, sampleRate, seed, samplesPerSymbol, signalPower, coefs, 2, 4, new double[] { 0, Math.PI / 4, Math.PI / 2, 3 * Math.PI / 4 }, new byte[] { 0, 1, 2, 3 });

			string filename = "RunFIRFilterSimulationObservableTest";
			filename = _testFileDumpDirec + filename + ".csv";
			var obs = mpsk.Subscribe(new MPSK_Modulation.IntermediateSimResultsReporter(filename));
			mpsk.RunFIRFilterSimulationObservable(500, 1700);
			obs.Dispose();
		}


		private void writeToCSV<T>(List<T> list, string filename, double sampleRate=0)
		{

			filename = _testFileDumpDirec + filename + ".csv";
			var csv = new StringBuilder();

			if (sampleRate != 0)
				csv.AppendLine(sampleRate.ToString());

			foreach (var item in list)
			{
				csv.AppendLine(item.ToString() + ",");
			}
			csv.Remove(csv.Length - 1, 1);

			File.WriteAllText(filename, csv.ToString());
			
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