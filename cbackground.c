#include <Python.h>
#include <import.h>
#include <stdio.h>

int main(
	int argc,
	char *argv[]
	//pyObject *py_func_Background_HubbleFunc
	)
{
	//initializations

	//Python objects (pointers):
	PyObject *py_module_Background, *py_class_Background, *py_tuple_Background_InitParams, *py_func_Background_Hubble;
	PyObject *py_instance_BackgroundClass, *py_dict_Background_ModuleAttributes, *py_float_Background_HubbleRate;

	//c objects
	double cosmo[5] = {.7, .12, 1.0e-5, .2, 2.0};
	Py_ssize_t ii;
	int check;


	// BEGIN CODE //
	printf("initialized variables\n");



	//begin Python wrapper
	Py_Initialize();
	printf("started python wrapper\n");

	//Load up the background module
	PySys_SetPath(".");
	py_module_Background = PyImport_ImportModuleNoBlock("background");
	printf("loaded python module\n");

	//Prepare the cosmology class init param tuple 
	py_tuple_Background_InitParams = PyTuple_New(5);
	printf("intialized python tuple\n");
	for (ii=0; ii < 5; ii++)
	{
		printf("putting %f at position %ld in tuple\n", cosmo[ii], ii);

		PyObject *py_temp_object = PyFloat_FromDouble(cosmo[ii]);

		check = PyTuple_SetItem(py_tuple_Background_InitParams, ii, py_temp_object);
		printf("PyTuple_SetItem returned %d\n", check);

		py_temp_object = PyTuple_GetItem(py_tuple_Background_InitParams, 0);
		if (!PyFloat_Check(py_temp_object)){printf("Tuple member isnt a float.  Why is that?\n");}
		else {printf("Found %f at position %ld in tuple\n", PyFloat_AsDouble(py_temp_object), ii);}
		Py_DECREF(py_temp_object);
	}
	printf("set python tuple\n");



	//Instance and initialize background module
	py_dict_Background_ModuleAttributes = PyModule_GetDict(py_module_Background);
	printf("grabbed module attributes\n");
	py_class_Background = PyDict_GetItemString(py_dict_Background_ModuleAttributes, "background");
	printf("grabbed background class\n");
	if (PyCallable_Check(py_class_Background))
	{
		py_instance_BackgroundClass = PyObject_CallObject(py_class_Background, NULL);
		if (!py_instance_BackgroundClass){printf("Bollocks.\n");}
	}
	else
	{
		return 1;
	}
	printf("instanced background class\n");
	PyObject *temptuple = PyTuple_New(1);
	check = PyTuple_SetItem(temptuple, 0,  PyFloat_FromDouble(2.0));
	printf("Created Tuple Object\n");
	if (PyObject_HasAttrString(py_instance_BackgroundClass, "hubble") == 1)
	{
		py_func_Background_Hubble = PyObject_GetAttrString(py_instance_BackgroundClass, "hubble");

		if (PyCallable_Check(py_func_Background_Hubble))
		{
			py_float_Background_HubbleRate = PyObject_Call(py_func_Background_Hubble, PyFloat_FromDouble(2.0), NULL);
			printf("called hubble method\n");
		}
		else{printf("nope, not a callable thing\n");}
	}
	if (py_float_Background_HubbleRate == NULL){printf("asfkl jqlwekurh !\n");}
	//printf("IS FLOAT returned %d", PyFloat_Check(py_float_Background_HubbleRate));

	// if (PyCallable_Check(py_func_Background_HubbleFunc))
	// {
	// 	PyObject *py_HubbleRate = PyObject_CallObject(py_func_Background_HubbleFunc, PyFloat_FromDouble(5));
	// 	printf("returned hubble rate of %f\n", PyFloat_AsDouble(py_HubbleRate));
	// }
	// else{printf("hubble_func not callable \n");}

	//Cleanup
	Py_DECREF(py_module_Background);
	Py_DECREF(py_class_Background);
	Py_DECREF(py_tuple_Background_InitParams);
	printf("DEALLOCATED MEMORY\n");

	//end Python wrapper
	Py_Finalize();
	printf("ended python wrapper\n");

	return 0;
}