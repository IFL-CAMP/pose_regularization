
#include <pose_regularization.h>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <matio.h>
#include <iostream>

using namespace std;

struct Settings {
    pose_regularization::Regularization p, q, r;
    double alpha, beta, inner_factor;
    int steps, inner_steps;
    string preset;
};

class MatFile {
public:
    explicit MatFile(const string& filename) {
        mat = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
        if (!mat) {
            throw std::runtime_error("Could not open file "+filename);
        }
    }

    vector<pose_regularization::SE3_ELEMENT> readPoseStreamMatrix(const string& varname) {

        matvar_t *matvar = Mat_VarReadInfo(mat, varname.c_str());
        if (!matvar) {
            throw std::runtime_error("Could not read variable "+varname);
        }

        int read_err = Mat_VarReadDataAll(mat, matvar);
        if (read_err) {
            throw std::runtime_error("Error reading data for "+varname);
        }
        if (matvar->rank != 3) {
            throw std::runtime_error("Unexpected rank "+std::to_string(matvar->rank));
        }
        if (matvar->data_type != MAT_T_DOUBLE) {
            throw std::runtime_error("Unexpected type "+std::to_string(matvar->data_type));
        }
        if (matvar->dims[0] != 4 || matvar->dims[1] != 4) {
            throw std::runtime_error("Unexpected row length "+std::to_string(matvar->dims[0]));
        }
//    std::cout << "Reading " << matvar->dims[2] << " rows" << std::endl;
        auto data_typed_ptr = reinterpret_cast<double*>(matvar->data);

        std::vector<pose_regularization::SE3_ELEMENT> ret(matvar->dims[2]);
        for (int i = 0; i < matvar->dims[2]; i++) {
            memcpy(ret[i].data(), data_typed_ptr+(i*16), 16*sizeof(double));
        }

//    Mat_VarPrint(field,1);

        Mat_VarFree(matvar);
        return ret;
    }

    int readIntVariable(const string& varname) {
        matvar_t *matvar = Mat_VarReadInfo(mat, varname.c_str());
        if (!matvar) {
            throw std::runtime_error("Could not read variable "+varname);
        }

        int read_err = Mat_VarReadDataAll(mat, matvar);
        if (read_err) {
            throw std::runtime_error("Error reading data for "+varname);
        }
        if (matvar->data_type != MAT_T_INT64) {
            throw std::runtime_error("Unexpected type "+std::to_string(matvar->data_type));
        }

        return reinterpret_cast<int*>(matvar->data)[0];
    }

    double readDoubleVariable(const string& varname) {
        matvar_t *matvar = Mat_VarReadInfo(mat, varname.c_str());
        if (!matvar) {
            throw std::runtime_error("Could not read variable "+varname);
        }

        int read_err = Mat_VarReadDataAll(mat, matvar);
        if (read_err) {
            throw std::runtime_error("Error reading data for "+varname);
        }
        if (matvar->data_type != MAT_T_DOUBLE) {
            throw std::runtime_error("Unexpected type "+std::to_string(matvar->data_type));
        }

        return reinterpret_cast<double*>(matvar->data)[0];
    }

    string readStringVariable(const string& varname) {
        matvar_t *matvar = Mat_VarReadInfo(mat, varname.c_str());
        if (!matvar) {
            throw std::runtime_error("Could not read variable "+varname);
        }

        int read_err = Mat_VarReadDataAll(mat, matvar);
        if (read_err) {
            throw std::runtime_error("Error reading data for "+varname);
        }
        if (matvar->class_type != MAT_C_CHAR) {
            throw std::runtime_error("Unexpected type "+std::to_string(matvar->data_type));
        }

        auto typed_data_ptr = reinterpret_cast<char*>(matvar->data);
        string ret(typed_data_ptr, typed_data_ptr+matvar->dims[1]);

        return ret;
    }

    void printVariableNames() {
        static char *mxclass[16] = {"cell", "struct", "object","char","sparse",
                                    "double","single","int8", "uint8","int16","uint16",
                                    "int32","uint32","int64","uint64","function"
        };
        matvar_t *matvar;
        size_t    nbytes;
        int       i;
        char size[32] = {'\0',};
        while ( NULL != (matvar = Mat_VarReadNextInfo(mat)) ) {
            printf("%-20s", matvar->name);
            if ( matvar->rank > 0 ) {
                int cnt = 0;
                printf("%8d", matvar->dims[0]);
                for ( i = 1; i < matvar->rank; i++ ) {
                    if ( ceil(log10(matvar->dims[i]))+1 < 32 )
                        cnt += sprintf(size+cnt,"x%d", matvar->dims[i]);
                }
                printf("%-10s",size);
            } else {
                printf("                    ");
            }
            printf("  %-18s\n",mxclass[matvar->class_type-1]);

            Mat_VarFree(matvar);
        }
    }

    Settings readSettings() {
        Settings ret {
            pose_regularization::Regularization(readIntVariable("p")),
            pose_regularization::Regularization(readIntVariable("q")),
            pose_regularization::Regularization(readIntVariable("r")),
            readDoubleVariable("alpha"),
            readDoubleVariable("beta"),
            readDoubleVariable("inner_factor"),
            readIntVariable("steps"),
            readIntVariable("inner_steps"),
            readStringVariable("preset"),
        };

        return ret;
    }

    ~MatFile() {
        Mat_Close(mat);
    }
private:
    mat_t *mat;
};

void repeatDenoising(const string& filename) {
    MatFile m{filename};

    vector<pose_regularization::SE3_ELEMENT> input = m.readPoseStreamMatrix("Poses_Noisy");
    vector<pose_regularization::SE3_ELEMENT> output(input.size());
    vector<pose_regularization::SE3_ELEMENT> expectedOutput = m.readPoseStreamMatrix("Poses_Regularized");
    Settings s = m.readSettings();

    pose_regularization::proximalPointAlgorithm(
            input, output,
            s.p, s.q, s.r,
            s.inner_factor,
            s.alpha, s.beta,
            s.steps, s.inner_steps
    );
    double ferr = 1e-5;
    for (int i = 0; i < expectedOutput.size(); i++)
        EXPECT_THAT(output[i], testing::Pointwise(testing::FloatNear(ferr), expectedOutput[i]));
}

TEST(Denoising, OpticalPreset)
{
    repeatDenoising("data/case5_preset_optical.mat");
}

TEST(Denoising, Em1Preset)
{
    repeatDenoising("data/case5_preset_em1.mat");
}

TEST(Denoising, Em2Preset)
{
    repeatDenoising("data/case5_preset_em2.mat");
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}