#include <Python.h>

int main(
	int argc,
	char *argv[]
	)

{
	PyObject *py_module_BackgroundLinker, *py_dict_BackgroundLinkerAttributes, *py_func_GetBackgroundHubble, *py_tuple_Function_Args;
	PyObject *py_func_Hubble, *py_float_HubbleAtZ, *py_tuple_Hubble_Args;
	int check;
	Py_ssize_t ii;
	double HubbleAtZ;

	Py_Initialize();
	//include local path
	PySys_SetPath(".");

	//import background_linker as py_module_BackgroundLinker
	py_module_BackgroundLinker = PyImport_ImportModuleNoBlock("background_linker");
	printf("imported background_linker\n");

	py_dict_BackgroundLinkerAttributes = PyModule_GetDict(py_module_BackgroundLinker);
	if (!PyDict_Contains(py_dict_BackgroundLinkerAttributes, PyString_FromString("get_background_hubble")))
	{
		printf("function 'get_background_hubble' not found \n");
		return 1;		
	}
	//py_funce_GetBackgroundHubble = py_module_BackgroundLinker.get_background_hubble
	py_func_GetBackgroundHubble = PyDict_GetItemString(py_dict_BackgroundLinkerAttributes, "get_background_hubble");
	printf("loaded get_background_hubble method\n");



	py_tuple_Function_Args = PyTuple_New(5);

	for (ii = 0; ii < 5; ii++)
	{
		PyObject *pValue = PyFloat_FromDouble(atof(argv[ii+1]));
		if (!pValue)
		{
			PyErr_Print();
			return 1;
		}
		check = PyTuple_SetItem(py_tuple_Function_Args, ii, pValue);
		if (check)
		{
			printf("Setting Tuple Value failed for location %ld\n", ii);
			return 1;
		}

	}
	printf("prepared function arguments\n");

	if (!PyCallable_Check(py_func_GetBackgroundHubble))
	{
		printf("'get_background_hubble' not callable\n");
		return 1;
	}
	
	//py_func_hubble = get_background_hubble(*py_tuple_function_args)
	py_func_Hubble = PyObject_Call(py_func_GetBackgroundHubble, py_tuple_Function_Args, NULL);
	printf("Called get_background_hubble\n");

	if (!PyCallable_Check(py_func_Hubble))
	{
		printf("Returned hubble object which is not callable\n");
		return 1;
	}

	//call hubble on some redshift to see if we succeeeded.
	py_tuple_Hubble_Args = PyTuple_New(1);
	check = PyTuple_SetItem(py_tuple_Hubble_Args, 0, PyFloat_FromDouble(15.0));
	if (check)
	{
		printf("Setting Tuple Value failed\n");
		return 1;
	}
	printf("Prepared hubble argument\n");

	py_float_HubbleAtZ = PyObject_CallObject(py_func_Hubble, py_tuple_Hubble_Args);
	if (!py_float_HubbleAtZ)
	{
		printf("Call to hubble function failed");
		return 1;
	}

	printf("Called hubble(15)\n");

	HubbleAtZ = PyFloat_AsDouble(py_float_HubbleAtZ);
	printf("We successully returned hubble(15) = %f\n", HubbleAtZ);
	Py_Finalize();
	return 0;
}
