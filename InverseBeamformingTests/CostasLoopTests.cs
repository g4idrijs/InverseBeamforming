using Microsoft.VisualStudio.TestTools.UnitTesting;
using InverseBeamforming;
using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace InverseBeamforming.Tests
{
	[TestClass()]
	public class CostasLoopTests
	{
		private string _testFileDumpDirec = Environment.GetEnvironmentVariable("USERPROFILE") + Path.DirectorySeparatorChar + @"Documents\TestFileDump\CostasLoopTests" + Path.DirectorySeparatorChar;

		[TestMethod()]
		public void runCostasLoopTest_NoGains()
		{
			double natfreq = 0;
			int order = 2;

			var cl = new CostasLoop(natfreq, order, false);

			Complex[] data = new Complex[100];

			for (int i = 0; i < data.Length; i++)
			{
				data[i] = 100 * new Complex(1, 0);
			}

			var output = cl.runCostasLoop(data);

			CollectionAssert.AreEqual(data, output);
		}

		[TestMethod()]
		public void runCostasLoopTest_PerfectData()
		{
			double natfreq = 0.25;
			int order = 2;
			Random rng = new Random();

			var cl = new CostasLoop(natfreq, order, false);

			Complex[] data = new Complex[100];

			for (int i = 0; i < data.Length; i++)
			{
				data[i] = new Complex(2*rng.Next(0,1)-1, 0);
			}

			var output = cl.runCostasLoop(data);

			CollectionAssert.AreEqual(data, output);
		}

		[TestMethod()]
		public void runCostasLoopTest_BPSKConvergenceTestWithStaticRotation()
		{
			double delta = .09;
			double natfreq = 0.25;
			int order = 2;
			Random rng = new Random();

			var cl = new CostasLoop(natfreq, order, false);

			Complex[] data = new Complex[100];
			Complex[] expected = new Complex[100];
			Complex rot = Complex.Exp(0.2 * Complex.ImaginaryOne);
			int N = 40;
			int num;
			for (int i = 0; i < data.Length; i++)
			{
				num = 2 * rng.Next(0, 2) - 1;
				expected[i] = new Complex(num, 0);
				data[i] = expected[i] * rot;
			}

			var output = cl.runCostasLoop(data);
			for (int i = N; i < expected.Length; i++)
			{
				Assert.AreEqual(expected[i].Magnitude, output[i].Magnitude,delta);
			}
			writeToCSV(expected, "BPSKConvergence_StaticRotation_Expected");
			writeToCSV(output, "BPSKConvergence_StaticRotation_output");
		}

		[TestMethod()]
		public void runCostasLoopTest_QPSKConvergenceTestWithStaticRotation()
		{
			double delta = .09;
			double natfreq = 0.25;
			int order = 4;
			Random rng = new Random();

			var cl = new CostasLoop(natfreq, order, false);

			Complex[] data = new Complex[100];
			Complex[] expected = new Complex[100];
			Complex rot = Complex.Exp(0.2 * Complex.ImaginaryOne);
			int N = 40;
			int num1,num2;
			for (int i = 0; i < data.Length; i++)
			{
				num1 = 2 * rng.Next(0, 2) - 1;
				num2 = 2 * rng.Next(0, 2) - 1;
				expected[i] = new Complex(num1, num2);
				data[i] = expected[i] * rot;
			}

			var output = cl.runCostasLoop(data);
			for (int i = N; i < expected.Length; i++)
			{
				Assert.AreEqual(expected[i].Magnitude, output[i].Magnitude, delta);
			}
			writeToCSV(expected, "QPSKConvergence_StaticRotation_Expected");
			writeToCSV(output, "QPSKConvergence_StaticRotation_output");
		}

		[TestMethod()]
		public void runCostasLoopTest_8PSKConvergenceTestWithStaticRotation()
		{
			double delta = .09;
			double natfreq = 0.25;
			int order = 8;
			Random rng = new Random();

			var cl = new CostasLoop(natfreq, order, false);

			Complex[] data = new Complex[100];
			Complex[] expected = new Complex[100];
			Complex[] constellation = new Complex[] { 0, Complex.Exp(Complex.ImaginaryOne*Math.PI / 4), Complex.Exp(Complex.ImaginaryOne * Math.PI / 2), Complex.Exp(Complex.ImaginaryOne * 3 * Math.PI / 4), Complex.Exp(Complex.ImaginaryOne * Math.PI), Complex.Exp(Complex.ImaginaryOne * 5 * Math.PI / 4), Complex.Exp(Complex.ImaginaryOne * 3 * Math.PI / 2), Complex.Exp(Complex.ImaginaryOne * 7 * Math.PI / 4 )};
			Complex rot = Complex.Exp(-Math.PI /(8 * Complex.ImaginaryOne));
			int N = 40;
			int num;
			for (int i = 0; i < data.Length; i++)
			{
				num = rng.Next(0, order);
				expected[i] = 2 * rot * constellation[num];
				data[i] = expected[i] * rot;
			}
			rot = Complex.Exp(0.1 * Complex.ImaginaryOne);
			var output = cl.runCostasLoop(data);
			for (int i = N; i < expected.Length; i++)
			{
				Assert.AreEqual(expected[i].Magnitude, output[i].Magnitude, delta);
			}
			writeToCSV(expected, "8PSKConvergence_StaticRotation_Expected");
			writeToCSV(output, "8PSKConvergence_StaticRotation_output");
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