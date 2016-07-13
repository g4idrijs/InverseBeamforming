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
	public class FilterTests
	{
		private string _testFileDumpDirec = Environment.GetEnvironmentVariable("USERPROFILE") + Path.DirectorySeparatorChar + @"Documents\TestFileDump\Filters" + Path.DirectorySeparatorChar;

		[TestMethod()]
		public void ExampleIIRCallTest()
		{
			string coeffFilename = _testFileDumpDirec + "IIR Filter.txt";
			string IIRFilename = _testFileDumpDirec + "IIR Filter Output.txt";
			string FFTFilename = _testFileDumpDirec + "FFT Output.txt";
			Filter.ExampleIIRCall(coeffFilename, FFTFilename, IIRFilename);
		}

		[TestMethod()]
		public void OtherFilterTest()
		{
			double[] a, b, q, u, v, w, r, z, input, outputVector;


			u = new double[] { 1, 1, 1 };
			v = new double[] { 1, 1, 0, 0, 0, 1, 1 };
			w = Filter.OtherFilter.conv(u, v);
			writeToCSV(w, "conv");

			a = new double[] { 1, 2, 3, 4 };
			b = new double[] { 10, 40, 100, 160, 170, 120 };
			q = Filter.OtherFilter.deconv(b, a);
			writeToCSV(q, "deconv");

			r = Filter.OtherFilter.deconvRes(b, a);
			writeToCSV(r, "deconvRes");

			a = new double[] { 2, -2.5, 1 };
			b = new double[] { 0.1, 0.1 };
			u = new double[31];
			for (int i = 1; i < u.Length; i++)
			{
				u[i] = 0.0;
			}
			u[0] = 1.0;
			z = Filter.OtherFilter.filter(b, a, u);
			writeToCSV(z, "filterImpulseResponse");

			a = new double[] { 1.0000, -3.518576748255174, 4.687508888099475, -2.809828793526308, 0.641351538057564 };
			b = new double[] { 0.020083365564211, 0, -0.040166731128422, 0, 0.020083365564211 };
			input = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

			outputVector = Filter.OtherFilter.filter(b, a, input);
			writeToCSV(outputVector, "filter");
		}

		private void writeToCSV<T>(List<T> list, string filename, double sampleRate = 0)
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