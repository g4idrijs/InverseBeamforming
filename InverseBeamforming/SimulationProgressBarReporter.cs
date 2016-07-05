using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	/// <summary>
	/// Holds methods that create and report on simulations of Bit error rate
	/// </summary>
	public partial class Simulations
	{
		/// <summary>
		/// Class that describes the process of being subscribed to a simulation
		/// </summary>
		public class SimulationProgressBarReporter : SimulationReporter
		{
			/// <summary>
			/// Event that indicates that the progress of the simulation has been updated
			/// </summary>
			public event SimulationProgressUpdated SimulationProgressBarUpdatedEvent;
			public delegate void SimulationProgressUpdated(SimulationProgressBarReporter s, ProgressUpdatedEventArgs e);

			/// <summary>
			/// Event that indicates that an error has occured in the simulation
			/// </summary>
			public event SimulationErrorOccured SimulationErrorOccuredEvent;
			public delegate void SimulationErrorOccured(SimulationProgressBarReporter s, ErrorOccuredEventArgs e);

			/// <summary>
			/// Percentage of the simulation that is completed
			/// </summary>
			public double Progress
			{
				get { return this._progress; }
			}
			protected double _progress;

			/// <summary>
			/// Length of time the simulation has been running, or if it is done, the length of time it was running.
			/// </summary>
			public TimeSpan LengthOfSimulation { get { return this._lengthOfSimulation; } }
			protected TimeSpan _lengthOfSimulation;

			/// <summary>
			/// Instance name
			/// </summary>
			protected string _instName;
			
			/// <summary>
			/// Construct a new instances of the SimulationReporter class
			/// </summary>
			/// <param name="simName">Name of the logfile to use for this reporter</param>
			public SimulationProgressBarReporter(string simName)
			{
				this._instName = simName;
				this._lengthOfSimulation = new TimeSpan(0, 0, 0);
			}
			
			/// <summary>
			/// Sets up the file to log the simulation
			/// </summary>
			public override void OnStart()
			{
				//Get the starting time
				startTime = DateTime.Now;

				//Reset the progress bar to 0% complete
				this.setProgressBar(0);

				//Tell subscribers that the progress bar was updated
				if(SimulationProgressBarUpdatedEvent!=null)
					SimulationProgressBarUpdatedEvent(this, new ProgressUpdatedEventArgs(0,new TimeSpan(0)));
			}

			/// <summary>
			/// Report the status of the simulation to the logfile
			/// </summary>
			/// <param name="sim">Simulation to report the status of</param>
			public override void OnNext(SingleSimulation sim)
			{
				//Get the intermediate results from the simulation
				IntermediateSimResults isr = sim.isr;

				//Update the progress bar
				this.setProgressBar(isr.PercentErrorHad);
				this._lengthOfSimulation = DateTime.Now - startTime;

				//Tell subscribers that the progress bar was updated
				if (SimulationProgressBarUpdatedEvent != null)
					SimulationProgressBarUpdatedEvent(this, new ProgressUpdatedEventArgs(isr.PercentErrorHad, this._lengthOfSimulation));
			}

			/// <summary>
			/// Report an error in the simulation to the logfile
			/// </summary>
			/// <param name="error">Exception describing the error</param>
			public override void OnError(Exception error)
			{
				//Tell subscribers that an error occured
				if (SimulationErrorOccuredEvent != null)
					SimulationErrorOccuredEvent(this, new ErrorOccuredEventArgs(error));
			}

			/// <summary>
			/// Write completion information to the logfile and unsuscribe from the simulation
			/// </summary>
			public override void OnCompleted()
			{
				//Get the end time of the simulation
				endTime = DateTime.Now;

				this.setProgressBar(100);
				this._lengthOfSimulation = endTime - startTime;

				//Tell subscribers that the progress bar was updated
				if (SimulationProgressBarUpdatedEvent != null)
					SimulationProgressBarUpdatedEvent(this, new ProgressUpdatedEventArgs(100,this._lengthOfSimulation));

				//Unsubscribe from the simulation
				this.Unsubscribe();
			}
			
			/// <summary>
			/// Private setter for the progress bar
			/// </summary>
			/// <param name="value">New value for the bar</param>
			private void setProgressBar(double value)
			{
				if (value >= 0 && value <= 100)
					this._progress = value;
				else
					throw new ArgumentOutOfRangeException("Progress", "Progress must be a percentage between 0 and 100 inclusive.");
			}
		}
	}

	/// <summary>
	/// Provides the basic progress of the simulation in an Event args inherited class
	/// </summary>
	public class ProgressUpdatedEventArgs : EventArgs
	{
		public double ProgressBar;
		public TimeSpan LengthOfSimulation;

		/// <summary>
		/// Constructs event args that hold the new progress of the simulation, and the time spent simulating so far
		/// </summary>
		/// <param name="percentageWorkDone">Percentage of the simulation that is done so far</param>
		/// <param name="LengthOfSimulationSoFar">Time spent simulating so far</param>
		public ProgressUpdatedEventArgs(double percentageWorkDone, TimeSpan LengthOfSimulationSoFar)
		{
			ProgressBar = percentageWorkDone;
			LengthOfSimulation = LengthOfSimulationSoFar;
		}
	}

	/// <summary>
	/// Event args class that can pass an exception
	/// </summary>
	public class ErrorOccuredEventArgs : EventArgs
	{
		public Exception e;

		public ErrorOccuredEventArgs(Exception e) : base()
		{
			this.e = e;
		}
	}
}
