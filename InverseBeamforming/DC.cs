using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/// <summary>
/// Contains all the components necesary to perform simulations regaurding digital communications systems
/// and inverse beamforming
/// </summary>
namespace InverseBeamforming
{
	/// <summary>
	/// Provides useful debugging methods, mainly used in the immediate window
	/// </summary>
	public static class DC
	{
		/// <summary>
		/// Displays only the differences between two byte arrays
		/// </summary>
		/// <param name="i">Input array</param>
		/// <param name="o">Output array</param>
		/// <returns>String containing all of the differences between the two byte arrays</returns>
		public static string DispDifferences(byte[] i, byte[] o)
		{
			StringBuilder s = new StringBuilder();

			for (int k = 0; k < i.Length; k++)
			{
				if (i[k] != o[k])
					s.Append("In: " + i[k] + " Out: " + o[k] + Environment.NewLine);
			}
			return s.ToString();
		}
	}
}
