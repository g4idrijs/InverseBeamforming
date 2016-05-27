using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	/// <summary>
	/// Defines an interface to control the way spreading codes are generated and handled.
	/// </summary>
	interface I_SpreadingCode
	{
		/// <summary>
		/// Holds the spreading code matrix
		/// </summary>
		double CodeMatrix { get; set; }

		/// <summary>
		/// Gets a full waveform length spreading code, ready to multiply with the signal. 
		/// </summary>
		/// <param name="user">Index of the spreading code to use (zero based)</param>
		/// <param name="numSamples">Number of samples per symbol duration</param>
		/// <param name="numChips">Number of chips per symbol duration</param>
		/// <param name="numSymbols">Number of symbols</param>
		/// <returns>Array containing the spreading code ready to be mixed with the signal</returns>
		double GetSpreadingCode(int user, int numSamples, int numChips, int numSymbols);
	}
}
