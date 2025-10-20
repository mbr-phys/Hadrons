#ifndef Hadrons_MDistil_Base_hpp_
#define Hadrons_MDistil_Base_hpp_


#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/DistilMatrix.hpp>
#include "Utils.hpp"

#define BASEGROUP "DistilMesonField"



using namespace Grid;
using namespace Hadrons;

template <typename T, typename Tio>
class ContractionDistilMesonField
{
private:
    std::map<std::string, DistilMesonFieldMatrix<Tio>> mf_;    //core 
public:
    ContractionDistilMesonField(std::string filename, unsigned int nt, TimerArray &timer);
    void load(void);
    std::string getName(void), getFileName(void);
    std::vector<std::vector<unsigned int>> getAvailTimeSources(void);
public:
    const DistilMesonFieldMatrix<T> operator()(const unsigned int t, 
        const unsigned int T1, const unsigned int T2) const
    {
        std::string key = "/" + std::to_string(t) + "/" + std::to_string(T1) + "-" + std::to_string(T2);
        return mf_.at(key).template cast <T>();
    }
private:
    TimerArray & tAr;
    const unsigned int nt_;
    std::string filename_;
    DistilMatrixIo<Tio> io_distil_;
};

template <typename T, typename Tio>
ContractionDistilMesonField<T,Tio>::ContractionDistilMesonField(std::string filename, unsigned int nt, TimerArray &timer)
: filename_(filename), nt_(nt), io_distil_(filename, BASEGROUP, nt), tAr(timer)
{
    this->load();
}

template <typename T, typename Tio>
std::string ContractionDistilMesonField<T,Tio>::getFileName()
{
    return filename_;
}

template <typename T, typename Tio>
void ContractionDistilMesonField<T,Tio>::load(void)
{
    unsigned int nDT = nt_;
    if(mf_.empty())
    {
        //assuming exact distillation with all time sources, assuming rho-phi type
        double total_time = 0.;
        for(int t=0; t<nt_; t++)
        {   
            double timer = 0.;
            //for(int T2=0; T2 < nDT; T2++)
            //{
            int T2 = t;
                double watch;
                std::string dataset_name = std::to_string(t) + "-" + std::to_string(T2);
                std::string key = "/" + std::to_string(t) + "/" + dataset_name;
                mf_.emplace(key, DistilMesonFieldMatrix<Tio>() );
                io_distil_.load(mf_.at(key), t, dataset_name, &watch, nullptr);
                timer += watch;
            //}
            CLOCK() << "= Read timeslice " << t << " at "<<  io_distil_.getSize()*nDT/timer*1.0e6/1024/1024 << " MB/s" << std::endl;
            total_time += timer;
        }
        double block_size = io_distil_.getSize();
        double timeslice_size = block_size * nDT;
        double total_size = timeslice_size * nt_;
        CLOCK() << "== Size per dilution block : " << block_size/1024/1024 << " MB" << std::endl;
        CLOCK() << "== Total size read : " << total_size/1024/1024 << " MB" << std::endl;
        CLOCK() << "== Average read speed : " << total_size/total_time*1.0e6/1024/1024 << " MB/s" << std::endl;
    }
}


#endif
