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
	}
}