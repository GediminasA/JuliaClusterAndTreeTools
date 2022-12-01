
def MAIN_DIR = "/mnt/jenkins/st16SseqJuliaTools"
echo "RUNNING WITH PARAMS: ${params.WITH_PARAMS}"
def BRANCH = ""

if(params.WITH_PARAMS){
    BRANCH = params.BRANCH
}
else{
    BRANCH = env.GIT_BRANCH
}
BRANCH = BRANCH.replaceAll(/[^\w]/, "_")

pipeline {
    agent {
        node {
            label 'master'
            customWorkspace "${MAIN_DIR}/st16SseqJuliaTools_$BRANCH"
        }
    }
    environment {
        JULIA_DIR="${MAIN_DIR}/JULIA_st16SseqJuliaTools"
        PATH_JULIA="$JULIA_DIR/julia-1.5.3/bin"
        PATH="$PATH_JULIA:$PATH"
        JULIA_PKGDIR="/mnt/jenkins/JULIA_PKGS"
        JULIA_DEPOT_PATH="$JULIA_DIR/.julia"
    }
    stages {
        stage('Install Julia') {
            steps {
                sh """#!/usr/bin/env bash
                    mkdir $JULIA_DIR
                    wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.3-linux-x86_64.tar.gz
                    tar -xvzf julia-1.5.3-linux-x86_64.tar.gz -C $JULIA_DIR
                """
            }
            post {
                success {
                    sh 'rm julia-1.5.3-linux-x86_64.tar.gz'
                }
            }
        }

        stage('Build Julia') {
            steps {
                sh """#!/usr/bin/env bash
                    bash ./setup.sh
                    julia -e 'using Pkg; Pkg.add("LightXML")'
                """
            }
        }

        stage('TESTS') {
            steps {
                sh """#!/usr/bin/env bash
                    julia --code-coverage=cov.info --project=. test/runtests.jl
                    julia test/julia_cov.jl cov.info julia_cov.xml \$PWD \$(ls -d \$PWD/src/*.jl)
                """
            }
        }
        stage('SonarQube') {
            steps{
                withSonarQubeEnv('SonarQube') {
                    script {
                        def sonarScanner = tool name: 'SQscanner', type: 'hudson.plugins.sonar.SonarRunnerInstallation'
                        sh """#!/usr/bin/env bash
                            ${sonarScanner}/bin/sonar-scanner -DskipTests -Dcobertura.skip
                        """
                    }
                }
            }
        }
    }
    post {
        cleanup{
            sh """#!/usr/bin/env bash
                rm $JULIA_DIR -R
            """
        }
    }
}
