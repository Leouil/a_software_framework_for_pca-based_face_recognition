
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lda.h"
#include "file_access.h"
//test code
#include <fstream>
#include <iostream>
#include "document.h"
#include "time.h"
using namespace std;
clock_t start, finish;
//test code

#include "file_access.h"

void init_param(struct corpus *cps,struct est_param* param,int topic_num);

// ./lda topic_num sample_num data model_name

int main(int argc, char *argv[])
{
    if ( argc < 5 )
    {
        printf("usage: ./lda topic_num sample_num data model_name\n");
        return 0;
    }
    int topic_num = atoi(argv[1]);
    int sample_num = atoi(argv[2]); // how many samples we want, we use these sample's to calculate phi and theta respectively and then use their averge to represent theta and phi
    char * data = argv[3];
    string model_name = argv[4];
    int burn_in_num = 1000; //how many iterations are there until we are out of burn-in period
    // set beta = 0.1 and alpha = 50 / K
    double beta = 0.1;
    double alpha = (double)50 / (double)topic_num;
    int SAMPLE_LAG = 10;
    struct corpus *cps;
    struct est_param param;
    cps = read_corpus(data);
    init_param(cps,&param,topic_num);
    cout << "parameter initialized" << endl;
    //do gibbs_num times sampling and compute averge theta,phi
    // burn in period
    //for (int iter_time=0; iter_time<burn_in_num + sample_num * SAMPLE_LAG; iter_time++)
    int iter_time = 0;
    int sample_time = 0;
    int p_size = topic_num<2?2:topic_num;
    double *p = (double*)malloc(sizeof(double)*p_size); //multinomial sampling tempoary storage space;
    double vbeta = cps->num_terms * beta; // v*beta is used in funciton sampling and this value never change so I put it here
    while (1)
    {
        start = clock();
        cout << "start sampling all document" << endl;
        //foreach documents, apply gibbs sampling
        for (int m=0; m<cps->num_docs; m++)
        {
            int word_index = 0; // word_index indicates that a word is the (word_index)th word in the document.
            double s_talpha = param.nd_sum[m] - 1 + topic_num * alpha; //nd_num[m]+topic_num * alpha will not change per document, so I put it here
            for (int l=0; l<cps->docs[m].length; l++)
            {
                for (int c=0; c<cps->docs[m].words[l].count; c++)
                {
                    param.z[m][word_index] = sampling(m,word_index,cps->docs[m].words[l].id,topic_num,cps,&param,alpha,beta,p,s_talpha,vbeta);
                    word_index++;
                }
            }
        }
        
        finish = clock();
        double dur = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "total time : " << dur << endl;
        
        
        if ((iter_time >= burn_in_num) && (iter_time % SAMPLE_LAG == 0))
        {
            calcu_param(&param, cps,topic_num,alpha,beta);
            cout << "\t" << sample_time+1 <<"# phi,theta calulation completed" << endl;
            sample_time++;
        }
        cout << iter_time+1 <<"# iteration completed"<<endl;
        iter_time++;
        if (sample_time ==sample_num)
        {
            break;
        }
    }
    // calculate average of phi and theta
    cout << "calculating the average of phi and theta" << endl;
    average_param(&param, cps,topic_num,alpha,beta,sample_num);
    cout << "parameter estimation completed, now saving model" << endl;
    save_model(cps,&param,model_name,alpha,beta,topic_num,sample_num);
    /*
     //test code
     ofstream out("ap.txt");
     
     for (int i=0; i<cps->num_docs; i++)
     {
     out << cps->docs[i].length << " ";
     for (int j=0; j < cps->docs[i].length; j++ )
     {
     out << cps->docs[i].words[j].id << ":" << cps->docs[i].words[j].count;
     if (j!=cps->docs[i].length-1)
     {
     out << " ";
     }
     }
     out << endl;
     }
     out.close();
     //test code
     */
    /*
     //test code
     for (int d=0; d<cps->num_docs; d++)
     {
     int word_index = 0; // word_index indicates that a word is the (word_index)th word in the document.
     int doc_length = 0;
     for (int l=0; l<cps->docs[d].length; l++)
     {
     
     for (int t=0; t<cps->docs[d].words[l].count; t++)
     {
     cout << param.z[d][word_index] << " ";
     word_index++;
     }
     }
     }
     //test code
     */
    return 0;
}

int sampling(int m, int n, int word_id, int topic_num, struct corpus* cps, struct est_param *param,double alpha, double beta,double* p,double s_talpha,double vbeta)
{
    int topic_id = param->z[m][n];
    param->nw[topic_id][word_id]--;
    param->nd[m][topic_id]--;
    param->nw_sum[topic_id]--;
    
    p[0] = (param->nw[0][word_id] + beta) / (param->nw_sum[0] + vbeta) * (param->nd[m][0] + alpha) / s_talpha;
    p[1] = p[0];
    for (int k = 1; k < topic_num-1; k++)
    {
        p[k] += (param->nw[k][word_id] + beta) / (param->nw_sum[k] + vbeta) * (param->nd[m][k] + alpha) / s_talpha;
        p[k+1] = p[k];
    }
    p[topic_num-1] += (param->nw[topic_num-1][word_id] + beta) / (param->nw_sum[topic_num-1] + vbeta) * (param->nd[m][topic_num-1] + alpha) / s_talpha;
    // choose a sample from topic
    double rand = ((double)random() / RAND_MAX) * p[topic_num - 1];
    int sample_topic;
    for (sample_topic = 0; sample_topic < topic_num; sample_topic++)
    {
        if (p[sample_topic] > rand)
        {
            break;
        }
    }
    // add new topic to statistics var
    param->nw[sample_topic][word_id]++;
    param->nd[m][sample_topic]++;
    param->nw_sum[sample_topic]++;
    return sample_topic;
}

void average_param(struct est_param *param, struct corpus* cps,int topic_num,double alpha,double beta, int sample_num)
{
    // theta
    for (int m=0; m<cps->num_docs; m++)
    {
        for (int k=0; k<topic_num; k++)
        {
            param->theta[m][k] /= sample_num;
        }
    }
    // phi
    for (int k=0; k<topic_num; k++)
    {
        for (int v=0; v<cps->num_terms; v++)
        {
            param->phi[k][v] /= sample_num;
        }
    }
}


void calcu_param(struct est_param *param, struct corpus* cps,int topic_num,double alpha,double beta)
{
    // theta
    for (int m=0; m<cps->num_docs; m++)
    {
        for (int k=0; k<topic_num; k++)
        {
            param->theta[m][k] += (param->nd[m][k] + alpha) / (param->nd_sum[m] + topic_num * alpha);
        }
    }
    // phi
    for (int k=0; k<topic_num; k++)
    {
        for (int v=0; v<cps->num_terms; v++)
        {
            param->phi[k][v] += (param->nw[k][v] + beta) / (param->nw_sum[k] + cps->num_terms * beta);
        }
    }
}

void set_zero(int **a,int rows,int columns)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<columns; j++)
        {
            a[i][j] = 0;
        }
    }
}

void set_zero_double(double **a,int rows,int columns)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<columns; j++)
        {
            a[i][j] = 0;
        }
    }
}


void init_param(struct corpus *cps, struct est_param* param,int topic_num)
{
    //1. alloc z[m][n], z[m][n] stands for topic assigned to nth word in mth document
    param->z = (int**)malloc (sizeof(int*) * cps->num_docs);
    for (int i=0; i<cps->num_docs; i++)
    {
        param->z[i] = (int*)malloc (sizeof(int) * cps->docs[i].num_term);
    }
    
    //2. alloc theta[m][k], theta[m][k] stands for the topic mixture proportion for document m
    param->theta = (double**)malloc (sizeof(double*) * cps->num_docs);
    for (int j=0; j<cps->num_docs; j++)
    {
        param->theta[j] = (double*)malloc (sizeof(double) * topic_num);
    }
    set_zero_double(param->theta,cps->num_docs,topic_num);
    
    //3. alloc phi[k][v], phi[k][v] stands for the probability of vth word in vocabulary is assigned to topic k
    param->phi = (double**)malloc (sizeof(double*) * topic_num);
    for (int z=0; z<topic_num; z++)
    {
        param->phi[z] = (double*)malloc (sizeof(double) * cps->num_terms);
    }
    set_zero_double(param->phi,topic_num,cps->num_terms);
    
    //4. alloc nd[m][k], nd[m][k] stands for the number of kth topic assigned to mth document
    param->nd = (int**)malloc (sizeof(int*) * cps->num_docs);
    for (int m=0; m< cps->num_docs; m++)
    {
        param->nd[m] = (int*)malloc(sizeof(int) * topic_num);
    }
    set_zero(param->nd, cps->num_docs, topic_num);
    
    //5. alloc nw[k][t], nw[k][t] stands for the number of kth topic assigned to tth term
    param->nw = (int**)malloc (sizeof(int*) * topic_num);
    for (int n=0; n < topic_num; n++)
    {
        param->nw[n] = (int*)malloc(sizeof(int) * cps->num_terms);
    }
    set_zero(param->nw, topic_num, cps->num_terms);
    
    //6. alloc nd_sum[m], nd_sum[m] total number of word in mth document
    param->nd_sum = (int*)malloc(sizeof(int)*cps->num_docs);
    memset(param->nd_sum,0,sizeof(int)*cps->num_docs);
    
    //7. alloc nw_sum[k], nw_sum[k] total number of terms assigned to kth topic
    param->nw_sum = (int*)malloc(sizeof(int)*topic_num);
    memset(param->nw_sum,0,sizeof(int)*topic_num);
    
    //8. init nd[m][k],nw[k][t],nd_sum[m],nw_sum[k]
    srandom(time(0)); // set seed for random function
    for (int d=0; d<cps->num_docs; d++)
    {
        int word_index = 0; // word_index indicates that a word is the (word_index)th word in the document.
        for (int l=0; l<cps->docs[d].length; l++)
        {
            for (int t=0; t<cps->docs[d].words[l].count; t++)
            {
                param->z[d][word_index] = (int)(((double)random() / RAND_MAX) * topic_num); // set z randomly
                param->nw[param->z[d][word_index]][cps->docs[d].words[l].id] ++;
                param->nd[d][param->z[d][word_index]]++;
                param->nw_sum[param->z[d][word_index]]++;
                
                word_index++;
            }
        }
        param->nd_sum[d] += word_index;
    }
}
