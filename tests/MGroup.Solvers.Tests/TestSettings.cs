namespace MGroup.Solvers.Tests
{
	using System;
	using System.IO;
	using System.Reflection;
	using System.Text.Json;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.Solvers.Tests.TempUtilityClasses;

	using Xunit;

	public static class TestSettings
	{
		private static readonly IReadOnlyList<IEnvironmentChoice> environmentChoicesToTest;

		// Explicit static constructor to tell C# compiler not to mark type as beforefieldinit. Only required for laziness.
		static TestSettings()
		{
			environmentChoicesToTest = new IEnvironmentChoice[]
			{
				new SequentialEnvironmentChoice(), new TplEnvironmentChoice(),
			};

			EnvironmentsToTestAsTheoryData = new TheoryData<IEnvironmentChoice>();
			EnvironmentsToTestAsTheoryData.Add(new SequentialEnvironmentChoice());
			EnvironmentsToTestAsTheoryData.Add(new TplEnvironmentChoice());

			LibsToTest = NativeLibsToTest.CreateWithNone();
			ProvidersToTestAsTheoryData = new TheoryData<IImplementationProvider>();
			var providersToTest = new List<IImplementationProvider>();
			ProvidersToTest = providersToTest;
			var providerChoicesToTest = new List<IImplementationProviderChoice>();
			ProviderChoicesToTest = providerChoicesToTest;

			// Always test the managed sequential implementation
			ProvidersToTestAsTheoryData.Add(new ManagedSequentialImplementationProvider());
			providersToTest.Add(new ManagedSequentialImplementationProvider());
			providerChoicesToTest.Add(new ManagedSequentialProviderChoice());

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
						ProvidersToTestAsTheoryData.Add(new NativeWin64ImplementationProvider());
						providersToTest.Add(new NativeWin64ImplementationProvider());
						providerChoicesToTest.Add(new NativeWin64ProviderChoice());
					}
				}
			}
			catch (Exception ex)
			{
				// If reading native lib options fails for any reason, do nothing (use only managed providers).
			}
		}

		public static NativeLibsToTest LibsToTest { get; }

		public static TheoryData<IEnvironmentChoice> EnvironmentsToTestAsTheoryData { get; }

		public static IReadOnlyList<IImplementationProvider> ProvidersToTest { get; }

		public static IReadOnlyList<IImplementationProviderChoice> ProviderChoicesToTest { get; }

		public static TheoryData<IImplementationProvider> ProvidersToTestAsTheoryData { get; }

		public const string SkipMessage =
			"This native library is not set to be tested. You can set it in NativeLibsToTest.json. See TestSettings.cs for more.";

		public static void CombineTheoryDataWithAllProviders<T1>(
			TheoryData<T1, IImplementationProviderChoice> theoryData, T1 param1)
		{
			foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
			{
				theoryData.Add(param1, provider);
			}
		}

		public static void CombineTheoryDataWithAllProviders<T1, T2>(
			TheoryData<T1, T2, IImplementationProviderChoice> theoryData, T1 param1, T2 param2)
		{
			foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
			{
				theoryData.Add(param1, param2, provider);
			}
		}

		public static void CombineTheoryDataWithAllProviders<T1, T2, T3>(
			TheoryData<T1, T2, T3, IImplementationProviderChoice> theoryData, T1 param1, T2 param2, T3 param3)
		{
			foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
			{
				theoryData.Add(param1, param2, param3, provider);
			}
		}

		public static void CombineTheoryDataWithAllProviders<T1, T2, T3, T4>(
			TheoryData<T1, T2, T3, T4, IImplementationProviderChoice> theoryData, T1 param1, T2 param2, T3 param3, T4 param4)
		{
			foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
			{
				theoryData.Add(param1, param2, param3, param4, provider);
			}
		}

		public static void CombineTheoryDataWithAllProviders<T1, T2, T3, T4, T5>(
			TheoryData<T1, T2, T3, T4, T5, IImplementationProviderChoice> theoryData,
			T1 param1, T2 param2, T3 param3, T4 param4, T5 param5)
		{
			foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
			{
				theoryData.Add(param1, param2, param3, param4, param5, provider);
			}
		}

		public static void CombineTheoryDataWithAllProviders<T1, T2, T3, T4, T5, T6>(
			TheoryData<T1, T2, T3, T4, T5, T6, IImplementationProviderChoice> theoryData,
			T1 param1, T2 param2, T3 param3, T4 param4, T5 param5, T6 param6)
		{
			foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
			{
				theoryData.Add(param1, param2, param3, param4, param5, param6, provider);
			}
		}

		public static void CombineTheoryDataWithAllProvidersAndEnvironments(
			TheoryData<IEnvironmentChoice, IImplementationProviderChoice> theoryData)
		{
			foreach (IEnvironmentChoice environment in environmentChoicesToTest)
			{
				foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
				{
					theoryData.Add(environment, provider);
				}
			}
		}

		public static void CombineTheoryDataWithAllProvidersAndEnvironments<T1>(
			TheoryData<T1, IEnvironmentChoice, IImplementationProviderChoice> theoryData, T1 param1)
		{
			foreach (IEnvironmentChoice environment in environmentChoicesToTest)
			{
				foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
				{
					theoryData.Add(param1, environment, provider);
				}
			}
		}

		public static void CombineTheoryDataWithAllProvidersAndEnvironments<T1, T2>(
			TheoryData<T1, T2, IEnvironmentChoice, IImplementationProviderChoice> theoryData, T1 param1, T2 param2)
		{
			foreach (IEnvironmentChoice environment in environmentChoicesToTest)
			{
				foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
				{
					theoryData.Add(param1, param2, environment, provider);
				}
			}
		}

		public static void CombineTheoryDataWithAllProvidersAndEnvironments<T1, T2, T3>(
			TheoryData<T1, T2, T3, IEnvironmentChoice, IImplementationProviderChoice> theoryData,
			T1 param1, T2 param2, T3 param3)
		{
			foreach (IEnvironmentChoice environment in environmentChoicesToTest)
			{
				foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
				{
					theoryData.Add(param1, param2, param3, environment, provider);
				}
			}
		}

		public static void CombineTheoryDataWithAllProvidersAndEnvironments<T1, T2, T3, T4>(
			TheoryData<T1, T2, T3, T4, IEnvironmentChoice, IImplementationProviderChoice> theoryData,
			T1 param1, T2 param2, T3 param3, T4 param4)
		{
			foreach (IEnvironmentChoice environment in environmentChoicesToTest)
			{
				foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
				{
					theoryData.Add(param1, param2, param3, param4, environment, provider);
				}
			}
		}

		public static void CombineTheoryDataWithAllProvidersAndEnvironments<T1, T2, T3, T4, T5>(
			TheoryData<T1, T2, T3, T4, T5, IEnvironmentChoice, IImplementationProviderChoice> theoryData,
			T1 param1, T2 param2, T3 param3, T4 param4, T5 param5)
		{
			foreach (IEnvironmentChoice environment in environmentChoicesToTest)
			{
				foreach (IImplementationProviderChoice provider in ProviderChoicesToTest)
				{
					theoryData.Add(param1, param2, param3, param4, param5, environment, provider);
				}
			}
		}
	}
}
