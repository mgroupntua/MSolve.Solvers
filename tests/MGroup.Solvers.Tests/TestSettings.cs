namespace MGroup.Solvers.Tests
{
	using System;
	using System.IO;
	using System.Reflection;
	using System.Text.Json;

	using MGroup.LinearAlgebra;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.LinearAlgebra.Implementations.NativeWin64.MKL;
	using MGroup.LinearAlgebra.Triangulation;

	using Xunit;

	public static class TestSettings
	{
		// Explicit static constructor to tell C# compiler not to mark type as beforefieldinit. Only required for laziness.
		static TestSettings()
		{
			LibsToTest = NativeLibsToTest.CreateWithNone();
			ProvidersToTest = new TheoryData<IImplementationProvider>();
			ProvidersToTest.Add(new ManagedSequentialImplementationProvider());

			try
			{
				// Read from JSON
				string execDirectory = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);
				string jsonFile = Path.Combine(execDirectory, "NativeLibsToTest.json");
				string jsonText = File.ReadAllText(jsonFile);

				// Currently there is no reason to test MKL without SuiteSparse. SuiteSparse dlls depend on MKL dlls.
				var options = new JsonSerializerOptions { PropertyNameCaseInsensitive = true };
				NativeLibsToTest? deserialized = JsonSerializer.Deserialize<NativeLibsToTest>(jsonText, options);
				if (deserialized != null)
				{
					LibsToTest = deserialized;
					if (LibsToTest.Win64IntelMkl && LibsToTest.Win64SuiteSparse)
					{
						ProvidersToTest.Add(new NativeWin64ImplementationProvider());
					}
				}
			}
			catch (Exception ex)
			{
				// If reading native lib options fails for any reason, do nothing (use only managed providers).
			}
		}

		public static NativeLibsToTest LibsToTest { get; }

		public static TheoryData<IImplementationProvider> ProvidersToTest { get; }

		public const string SkipMessage =
			"This native library is not set to be tested. You can set it in NativeLibsToTest.json. See TestSettings.cs for more.";
	}
}
