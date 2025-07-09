#ifndef HDF5_UTILS_H
#define HDF5_UTILS_H
#include <hdf5.h>
#include <mpi.h>
 
hid_t my_H5Gcreate(hid_t loc_id, const char *groupname, size_t size_hint);
herr_t my_H5Gclose(hid_t group_id, const char *groupname);
hid_t my_H5Fcreate(const char *fname, unsigned int flags, hid_t fcpl_id, hid_t fapl_id);
herr_t my_H5Fclose(hid_t file_id, const char *fname);
hid_t my_H5Fopen(const char *fname, unsigned int flags, hid_t fapl_id);
hid_t my_H5Screate(H5S_class_t type);
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
void write_vector_attribute(hid_t handle, const char *attr_name, const void *buf, hid_t mem_type_id, int length);
void write_scalar_attribute(hid_t handle, const char *attr_name, const void *buf, hid_t mem_type_id);
herr_t my_H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count,
                              const hsize_t *block);
herr_t my_H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buf,
                   const char *datasetname);
herr_t my_H5Dclose(hid_t dataset_id, const char *datasetname);
hid_t my_H5Screate_simple(int rank, const hsize_t *current_dims, const hsize_t *maximum_dims);
hid_t my_H5Dcreate(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id);

hid_t my_H5Pcreate(hid_t cls_id);
herr_t my_H5Pclose(hid_t plist_id);
herr_t my_H5Pset_fapl_mpio(hid_t plist_id, MPI_Comm comm, MPI_Info info);
herr_t my_H5Pset_dxpl_mpio(hid_t plist_id, H5FD_mpio_xfer_t xfer_mode);
#endif
